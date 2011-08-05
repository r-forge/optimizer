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
    f = expr2fun(par.name, expr);
    f2 = function(x) unlist(f(x));
    n = length(expr);
    result = character(n);
    test.box1 = sweep(diag(runif(n)), 2, par.val, "+");
    test.box2 = sweep(diag(runif(n)), 2, par.val, "+");
    test.box = rbind(test.box1, test.box2);
    val.box = apply(test.box, 1, f2);
    res.box = apply(val.box, 1, function(x) length(unique(x)));
    val.lin = apply(test.box + runif(1), 1, f2);
    res.lin = val.lin - val.box;
    test.lin = apply(res.lin, 1, function(x) length(unique(x)));
    result[test.lin < 1.5 & res.box < 3.5] = "box";
    result[test.lin < 1.5 & res.box > 3.5] = "linear";
    result[result == ""] = "nonlinear";
    return(result);
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
    }
    constr = gsub("^\\{|\\}$", "", src[2]);
    constr = tidyCode(constr);
    constr = tidyConstr(constr);
    withinBound = checkBound(constr$expr, parName, parVal);
    constrType = classifyConstr(constr$expr, parName, parVal);
    return(list(assignment = src[1], parName = parName, parVal = parVal,
                constrType = constrType, withinBound = withinBound));
}

# src = "x = c(0, 1, 2, 3);
#        
#        {
#            x[1] >= 0;
#            x[2] <= 3;
#            x[1] + x[2] <= 5 * x[3];
#            x[3]^2 + exp(x[4]) <= 10;
#        }";
# reportPar(src);
