tidyCode = function(src)
{
    exprs = parse(text = src);
    n = length(exprs);
    res = character(n);
    for (i in 1:n)
    {
        dep = paste(deparse(exprs[i]), collapse = "\n");
        res[i] = substring(dep, 12, nchar(dep) - 1);
    }
    return(res);
}

analyzeAssignment = function(src)
{
    env = new.env();
    eval(parse(text = src), envir = env);
    parName = ls(envir = env);
    if(!length(parName))
    {
        parName = "";
        parVal = "";
    }else{
        parVal = get(parName, envir = env);
    }
    return(list(parName = parName, parVal = parVal));
}

tidyConstr = function(src)
{
    exprs = gsub(" < ", " <= ", src);
    exprs = gsub(" > ", " >= ", exprs);
    le = grepl(" <= ", exprs);
    ge = grepl(" >= ", exprs);
    exprs = gsub("(<=|==|>=)", "# \\1 #", exprs);
    exprs = strsplit(exprs, " # ");
    newexpr = character(length(exprs));
    newexpr[le] = sapply(exprs[le], function(x) sprintf("%s - (%s)",
                                                        x[1], x[3]));
    newexpr[ge] = sapply(exprs[ge], function(x) sprintf("%s - (%s)",
                                                        x[3], x[1]));
    return(list(expr = parse(text = newexpr), src = paste(newexpr, "<= 0")));
}

expr2fun = function(par.name, expr)
{
    f = function(x)
    {
        envir = environment();
        assign(par.name, x, envir = envir);
        return(lapply(expr, eval, envir = envir));
    }
    return(f);
}

checkBound = function(expr, par.name, par.val)
{
    if(!length(expr)) return(NULL);
    f = expr2fun(par.name, expr);
    val = f(par.val);
    withinBound = sapply(val, function(x) all(x <= 0));
    return(withinBound);
    #index = (1:length(val))[!withinBound];
    #if(!length(index)) index = 0;
    #return(index);
}

classifyConstr = function(expr, par.name, par.val)
{
    if(!length(expr)) return(NULL);
    f = expr2fun(par.name, expr);
    f2 = function(x) unlist(f(x));
    n = length(expr);
    p = length(par.val);
    result = character(n);
    test.box1 = sweep(diag(runif(p)), 2, par.val, "+");
    test.box2 = sweep(diag(runif(p)), 2, par.val, "+");
    # Each row is a vector to be evaluated by f2
    test.box = rbind(test.box1, test.box2);
    val.box = apply(test.box, 1, f2);
    dim(val.box) = c(n, 2 * p);
    val.box = t(val.box);
    res.box = apply(val.box, 2, function(x) length(unique(x)));
    val.lin = apply(test.box + runif(1), 1, f2);
    dim(val.lin) = c(n, 2 * p);
    val.lin = t(val.lin) - val.box;
    res.lin = apply(val.lin, 2, function(x) length(unique(x)));
    result[res.lin < 1.5 & res.box < 3.5] = "box";
    result[res.lin < 1.5 & res.box > 3.5] = "linear";
    result[result == ""] = "nonlinear";
    
    varAssign = function(var, varName)
    {
        if(is.matrix(var))
        {
            assignExpr = paste(as.numeric(round(var, 10)), collapse = ",");
            assignExpr = sprintf("%s <- matrix(c(%s), %s, %s)", varName,
                                 assignExpr, nrow(var), ncol(var));
        } else {
            assignExpr = sprintf("%s <- %s", varName, paste(deparse(var),
                                                            collapse = "\n"));
        }
        assignExpr = tidyCode(assignExpr);
        return(assignExpr);
    }
    
    # Get the lower and upper bound of box constraint
    # lb <= x <= ub
    boxConstrIndex = (1:n)[result == "box"];
    if(!length(boxConstrIndex))
    {
        a = bound = lb = ub = boxConstrAssign = NULL;
    } else {
        boxElemIndex = apply(val.box[, boxConstrIndex, drop = FALSE], 2,
                             function(x) which(!duplicated(x, fromLast = TRUE))[1]);
        x1 = test.box[1, boxElemIndex];
        y1 = val.box[1, boxConstrIndex];
        x2 = test.box[2, boxElemIndex];
        y2 = val.box[2, boxConstrIndex];
        a = (y1 - y2) / (x1 - x2);
        bound = x1 - a * y1;
        lb = rep(-Inf, p);
        ub = rep(Inf, p);
        lb[boxElemIndex[a < 0]] = bound[a < 0];
        ub[boxElemIndex[a > 0]] = bound[a > 0];
        boxConstrAssign = c(varAssign(lb, "LB"), varAssign(ub, "UB"));
    }
    # Get the coefficient matrix and vector of linear constraint
    # A %*% x - b <= 0
    linConstrIndex = (1:n)[result == "linear"];
    if(!length(linConstrIndex))
    {
        A = b = linConstrAssign = NULL;
    } else {
        x = cbind(test.box[1:(p + 1), ], -1);
        m = solve(x, val.box[1:(p + 1), linConstrIndex]);
        dim(m) = c(1 + p, length(linConstrIndex));
        A = round(t(m[1:p, ]), 10);
        b = round(m[p + 1, ], 10);
        linConstrAssign = c(varAssign(A, "linMat"), varAssign(b, "linUB"));
    }
    # Combine box and linear constraints
    boxlinConstrIndex = c(boxConstrIndex, linConstrIndex);
    if(!length(boxlinConstrIndex))
    {
        A2 = b2 = boxlinConstrAssign = NULL;
    } else {
        x = cbind(test.box[1:(p + 1), ], -1);
        m = solve(x, val.box[1:(p + 1), boxlinConstrIndex]);
        dim(m) = c(1 + p, length(boxlinConstrIndex));
        A2 = round(t(m[1:p, ]), 10);
        b2 = round(m[p + 1, ], 10);
        boxlinConstrAssign = c(varAssign(A2, "linMat"), varAssign(b2, "linUB"));
    }
    # Get the function and upper bound of nonlinear constraint
    # f(x) <= fb
    nonlinConstrIndex = (1:n)[result == "nonlinear"];
    nnonlin = length(nonlinConstrIndex);
    if(!nnonlin)
    {
        ineqFun = ineqUB = nonlinearConstrAssign = NULL;
        ineqFun2 = ineqUB2 = nonlinearConstrAssign2 = NULL;
    } else {
        nonlinExprs = as.character(expr)[nonlinConstrIndex];
        assignExprs = paste("v", 1:nnonlin, "<-", nonlinExprs, sep = "", collapse = "\n");
        varNames = paste("v", 1:nnonlin, sep = "", collapse = ",");
        returnExprs = sprintf("return(c(%s))", varNames);
        ineqFunBody = sprintf("{%s \n %s}", assignExprs, returnExprs);
        ineqFun = function(param) NULL;
        names(formals(ineqFun)) = par.name;
        body(ineqFun) = parse(text = ineqFunBody);
        ineqUB = rep(0, nnonlin);
        nonlinearConstrAssign = c(varAssign(ineqFun, "ineqFun"),
                                  varAssign(ineqUB, "ineqUB"));
        # Combine linear and nonlinear constraints
        if(!is.null(linConstrAssign))
        {
            linExprs = sprintf("\n%s\n%s\n", linConstrAssign[1], linConstrAssign[2]);
            linConstrExprs = sprintf("linConstr <- linMat %s %s - linUB", "%*%", par.name);
            returnExprs2 = sprintf("return(c(linConstr, %s))", varNames);
            ineqFunBody2 = sprintf("{%s\n %s\n %s \n %s}", linExprs,
                                   linConstrExprs, assignExprs, returnExprs2);
            ineqFun2 = ineqFun;
            body(ineqFun2) = parse(text = ineqFunBody2);
            ineqUB2 = rep(0, nnonlin + length(linConstrIndex));
            nonlinearConstrAssign2 = c(varAssign(ineqFun2, "ineqFun"),
                                       varAssign(ineqUB2, "ineqUB"));
        } else {
            ineqFun2 = ineqFun;
            ineqUB2 = ineqUB;
            nonlinearConstrAssign2 = nonlinearConstrAssign;
        }
    }
    return(list(type = result,
                box = list(lb = lb, ub = ub, assign = boxConstrAssign),
                linear = list(A = A, b = b, assign = linConstrAssign),
                boxlinear = list(A = A2, b = b2, assign = boxlinConstrAssign),
                nonlinear = list(ineqFun = ineqFun, ineqUB = ineqUB,
                                 assign = nonlinearConstrAssign),
                nonlinear2 = list(ineqFun = ineqFun2, ineqUB = ineqUB2,
                                  assign = nonlinearConstrAssign2)));
}

reportPar = function(src)
{
    src = tidyCode(src);
    assignment = src[1];
    assignment = analyzeAssignment(assignment);
    parName = assignment$parName;
    parVal = assignment$parVal;
    if(length(src) < 2)
    {
        constrType = NULL;
        withinBound = NULL;
        boxConstr = NULL;
        linearConstr = NULL;
        boxlinearConstr = NULL;
        nonlinearConstr = NULL;
        nonlinearConstr2 = NULL;
    }else{
        constr = gsub("^.*\\{|\\}$", "", src[2]);
        constr = tidyCode(constr);
        constr = tidyConstr(constr);
        withinBound = checkBound(constr$expr, parName, parVal);
        constrClass = classifyConstr(constr$expr, parName, parVal);
        constrType = constrClass$type;
        boxConstr = constrClass$box;
        linearConstr = constrClass$linear;
        boxlinearConstr = constrClass$boxlinear;
        nonlinearConstr = constrClass$nonlinear;
        nonlinearConstr2 = constrClass$nonlinear2;
    }
    return(list(assignment = src[1], parName = parName, parVal = parVal,
                constrType = constrType, withinBound = withinBound,
                boxConstr = boxConstr, linearConstr = linearConstr,
                boxlinearConstr = boxlinearConstr,
                nonlinearConstr = nonlinearConstr,
                nonlinearConstr2 = nonlinearConstr2));
}

# src = "x = c(0, 1, 2, 3);
# 
# constr = {
#            x[2] >= -1;
#            x[1] <= 3;
#            x[2] + 1 <= 5;
#            x[1] + x[2] <= 5 * x[3];
#            x[2] - x[1] >= -1;
#            x[3]^2 + exp(x[4]) <= 10;
#            x[4]^2 >= 3;
#          }";
# 
# src = tidyCode(src);
# assignment = src[1];
# par.name = analyzeAssignment(assignment)$parName;
# par.val = analyzeAssignment(assignment)$parVal;
# constr = gsub("^.*\\{|\\}$", "", src[2]);
# constr = tidyCode(constr);
# constr = tidyConstr(constr);
# expr = constr$expr;
# reportPar(src);


