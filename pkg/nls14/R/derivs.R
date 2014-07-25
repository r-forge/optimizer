# R-based replacement for deriv() function
#
#  -- needs a database of derivative definitions, extendable by the user
#  -- implement this as an environment, with keys being the name of the function whose
#     derivative is being defined.
#  -- have newDeriv function that adds a new derivative, e.g.

# myDerivs <- new.env(parent = sysDerivs)
# newDeriv(log(x, base = exp(1)), ..., derivs = myDerivs)
# etc.

# Need simplify function that does what the C simplification does, but on a language object.

sysDerivs <- new.env(parent = emptyenv())

newDeriv <- function(expr, deriv, envir = new.env(parent=sysDerivs)) {
    expr <- substitute(expr)
    deriv <- substitute(deriv)
    fn <- as.character(expr[[1]])
    args <- expr[-1]
    argnames <- names(args)
    if (is.null(argnames))
    	argnames <- rep("", length(args))
    required <- which(argnames == "")
    for (i in required)
    	argnames[i] <- as.character(args[[i]])    	
    assign(fn, list(head=expr, argnames = argnames, 
                    required = required, deriv=deriv), envir = envir)
    invisible(envir[[fn]])
}
         
Deriv <- function(expr, name, envir = sysDerivs, do_substitute = TRUE) {
    Recurse <- function(expr) {
    	if (is.call(expr)) {
    	    if (as.character(expr[[1]]) == "D")
    	    	expr <- Deriv(expr[[2]], name, envir, do_substitute = FALSE)
    	    else
    	    	for (i in seq_along(expr)[-1])
    	    	    expr[[i]] <- Recurse(expr[[i]])
    	}
    	expr
    }
    Simplify <- function(expr) {
    	if (is.call(expr)) {
    	    for (i in seq_along(expr)[-1])
    	        expr[[i]] <- Simplify(expr[[i]])
    	    fn <- as.character(expr[[1]])
    	    switch(fn,
    	    "+" = if (length(expr) == 2 || 
    	              (is.numeric(expr[[3]]) && expr[[3]] == 0))
    	    	     expr[[2]]
    	    	  else if (is.numeric(expr[[2]]) && expr[[2]] == 0)
    	    	     expr[[3]]
    	    	  else
    	    	     expr,
    	    "-" = if (length(expr) == 2)
    	    	      expr
    	    	  else if (is.numeric(expr[[3]]) && expr[[3]] == 0)
    	    	     expr[[2]]
    	    	  else if (is.numeric(expr[[2]]) && expr[[2]] == 0)
    	    	     call("-", expr[[3]])
    	    	  else
    	    	     expr,     
    	    "*" = if (is.numeric(expr[[3]]) && expr[[3]] == 1)
    	    	     expr[[2]]
    	    	  else if (is.numeric(expr[[2]]) && expr[[2]] == 1)
    	    	     expr[[3]]
    	    	  else if (is.numeric(expr[[2]]) && expr[[2]] == 0)
    	    	     0    	    	     
    	    	  else if (is.numeric(expr[[3]]) && expr[[3]] == 0)
    	    	     0
    	    	  else
    	    	     expr,
    	    "/" = if (is.numeric(expr[[3]]) && expr[[3]] == 1)
    	    	     expr[[2]]
    	    	  else
    	    	     expr,
    	    "(" = expr[[2]],	     
    	    	expr)
    	} else
    	    expr
    }
    if (do_substitute)
    	expr <- substitute(expr)
    if (is.expression(expr))
    	return(as.expression(lapply(expr, Deriv, name=name)))
    else if (is.numeric(expr))
    	return(0)
    else if (is.name(expr))
    	if (as.character(expr) == name)
    	    return(1)
    	else
    	    return(0)
    else if (is.call(expr)) {
    	fn <- as.character(expr[[1]])

    	model <- envir[[fn]]
    	if (is.null(model))
    	    stop("no derivative known for '", fn, "'")
   
        args <- expr[-1]
        argnames <- names(args)
        if (is.null(argnames)) 
            argnames <- rep("", length(args))
        modelnames <- model$argnames
        argnum <- pmatch(argnames, modelnames)
        if (any(bad <- is.na(argnum) & argnames != ""))
            stop("Argument names not matched: ", paste(argnames[bad], collapse=", "))
        unused <- setdiff(seq_along(modelnames), argnum)
        nonamecount <- sum(is.na(argnum))
        length(unused) <- nonamecount
        argnum[which(is.na(argnum))] <- unused
        default <- setdiff(seq_along(modelnames), argnum)
        if (length(bad <- setdiff(model$required, argnum)))
            stop("Missing required arguments: ", paste(modelnames[bad], collapse=", "))
        
        # Now do the substitutions
        subst <- list()
        subst[argnum] <- as.list(args)
        subst[default] <- as.list(model$expr[-1])[default]
        names(subst) <- modelnames
        result <- do.call(substitute, list(model$deriv, subst))
        result <- Recurse(result)
        Simplify(result)
    }        
}

   
newDeriv(log(x), 
         (1/x)*D(x), envir = sysDerivs)
newDeriv(x + y, D(x) + D(y), envir = sysDerivs)
newDeriv(x*y, x*D(y) + D(x)*y, envir = sysDerivs)

#print(Deriv(log(y), "y"))
#print(Deriv(log(1+z), "z"))
#print(Deriv(log(x)*log(y), "x"))

    