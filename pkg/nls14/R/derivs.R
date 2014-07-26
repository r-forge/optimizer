# R-based replacement for deriv() function

sysDerivs <- new.env(parent = emptyenv())

newDeriv <- function(expr, deriv, envir = sysDerivs) {
    if (missing(expr))
    	return(ls(envir))
    expr <- substitute(expr)
    fn <- as.character(expr[[1]])
    if (missing(deriv)) 
    	return(envir[[fn]])
    deriv <- substitute(deriv)
    args <- expr[-1]
    argnames <- names(args)
    if (is.null(argnames))
    	argnames <- rep("", length(args))
    required <- which(argnames == "")
    for (i in required)
    	argnames[i] <- as.character(args[[i]])    	
    assign(fn, list(expr=expr, argnames = argnames, 
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

    if (do_substitute)
    	expr <- substitute(expr)
    if (is.expression(expr))
    	return(as.expression(lapply(expr, Deriv, name=name, do_substitute = FALSE)))
    else if (is.numeric(expr) || is.logical(expr))
    	return(0)
    else if (is.name(expr))
    	if (as.character(expr) == name)
    	    return(1)
    	else
    	    return(0)
    else if (is.call(expr)) {
    	fn <- as.character(expr[[1]])
	if (fn == "expression")
	    return(as.expression(lapply(as.list(expr)[-1], Deriv, name=name, do_substitute = FALSE)))
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

Simplify <- function(expr) {
    if (is.expression(expr))
    	return(as.expression(lapply(expr, Simplify)))    
    
    isFALSE <- function(x) identical(FALSE, x)
    isZERO <- function(x) is.numeric(x) && length(x) == 1 && x == 0
    isONE  <- function(x) is.numeric(x) && length(x) == 1 && x == 1
    isMINUSONE <- function(x) is.numeric(x) && length(x) == 1 && x == -1

    if (is.call(expr)) {
	for (i in seq_along(expr)[-1])
	    expr[[i]] <- Simplify(expr[[i]])
	fn <- as.character(expr[[1]])
	if (length(expr) == 2) 
	    switch(fn,
	    "+" = expr[[2]],
	    "-" = if (is.numeric(expr[[2]]))
		      -expr[[2]]
		  else if (is.call(expr[[2]]) && length(expr[[2]]) == 2 && as.character(expr[[c(2,1)]]) == "-")
		      expr[[c(2,2)]]
		  else 
		     expr,     
	    "exp" = if (is.call(expr[[2]]) && as.character(expr[[c(2,1)]]) == "log"
					   && length(expr[[2]]) == 2)
			expr[[c(2,2)]]
		    else
			expr,
	    "log" = if (is.call(expr[[2]]) && as.character(expr[[c(2,1)]]) == "exp")
			expr[[c(2,2)]]
		    else if (is.numeric(expr[[2]]))
		    	log(expr[[2]])
		    else
			expr,
	    "!" =   if (isTRUE(expr[[2]]))
	    		FALSE
	            else if (isFALSE(expr[[2]]))
	            	TRUE
	            else
	            	expr,
	    "(" = expr[[2]],	     
		expr)
	else if (length(expr) == 3)
	    switch(fn,
	    "+" = if (isZERO(expr[[3]]))
		     expr[[2]]
		  else if (isZERO(expr[[2]]))
		     expr[[3]]
		  else if (is.numeric(expr[[2]]) && is.numeric(expr[[3]]))
		     expr[[2]] + expr[[3]]
		  else
		     expr,
	    "-" = if (isZERO(expr[[3]]))
		     expr[[2]]
		  else if (isZERO(expr[[2]]))
		     Simplify(call("-", expr[[3]]))
		  else if (is.numeric(expr[[2]]) && is.numeric(expr[[3]]))
		     expr[[2]] - expr[[3]]
		  else
		     expr,     
	    "*" = if (isONE(expr[[3]]))
		     expr[[2]]
		  else if (isONE(expr[[2]]))
		     expr[[3]]
		  else if (isMINUSONE(expr[[3]]))
		     Simplify(call("-", expr[[2]]))
		  else if (isMINUSONE(expr[[2]]))
		     Simplify(call("-", expr[[3]]))
		  else if (isZERO(expr[[2]]))
		     0    	    	     
		  else if (isZERO(expr[[3]]))
		     0
		  else if (is.numeric(expr[[2]]) && is.numeric(expr[[3]]))
		     expr[[2]] * expr[[3]]
		  else
		     expr,
	    "/" = if (isONE(expr[[3]]))
		     expr[[2]]
		  else if (isMINUSONE(expr[[3]]))
		     Simplify(call("-", expr[[2]]))
		  else if (isZERO(expr[[2]]))
		     0
		  else if (is.numeric(expr[[2]]) && is.numeric(expr[[3]]))
		     expr[[2]] / expr[[3]]
		  else
		     expr,
	    "^" = if (isONE(expr[[3]]))
		     expr[[2]]
		  else if (is.numeric(expr[[2]]) && is.numeric(expr[[3]]))
		     expr[[2]] ^ expr[[3]]
		  else 
		     expr,
	    "log" = if (is.call(expr[[2]]) && as.character(expr[[c(2,1)]]) == "exp")
			Simplify(call("/", expr[[c(2,2)]], call("log", expr[[3]])))
		    else
			expr,	
	    "&&" = if ( isFALSE(expr[[2]]) || isFALSE(expr[[3]]))
	               FALSE
	           else if (isTRUE(expr[[2]]))
	               expr[[3]]
	           else if (isTRUE(expr[[3]]))
	               expr[[2]]
	           else
	               expr,
	    "||" = if (  isTRUE(expr[[2]]) || isTRUE(expr[[3]]))
	               TRUE
	           else if (isFALSE(expr[[2]]))
	               expr[[3]]
	           else if (isFALSE(expr[[3]]))
	               expr[[2]]
	           else
	               expr,
	    "if" = if (isTRUE(expr[[2]]))
	    	       expr[[3]]
	    	   else if (isFALSE(expr[[2]]))
	    	       NULL
	    	   else
	    	       expr,   
		expr)
	else if (length(expr) == 4)
	    switch(fn,
    	    "if" = if (isTRUE(expr[[2]]))
    	    	       expr[[3]]
    	    	   else if (isFALSE(expr[[2]]))
    	    	       expr[[4]]
    	    	   else if (identical(expr[[3]], expr[[4]]))
    	    	       expr[[3]]
    	    	   else
    	    	       expr,
    	    expr)
	else
	    expr 
    } else
	expr
}

# These are the derivatives supported by deriv()

newDeriv(log(x, base = exp(1)), 
         D(x)/(x*log(base)) - (1/base)*log(x, base)/log(base)*D(base), envir = sysDerivs)
newDeriv(exp(x), exp(x)*D(x), envir = sysDerivs)
newDeriv(sin(x), cos(x)*D(x), envir = sysDerivs)
newDeriv(cos(x), -sin(x)*D(x), envir = sysDerivs)
newDeriv(tan(x), 1/cos(x)^2*D(x), envir = sysDerivs)
newDeriv(sinh(x), cosh(x)*D(x), envir = sysDerivs)
newDeriv(cosh(x), sinh(x)*D(x), envir = sysDerivs)
newDeriv(sqrt(x), D(x)/2/sqrt(x), envir = sysDerivs)
newDeriv(pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),
  (if (lower.tail && !log.p) dnorm((q-mean)/sd)*(D(x)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (!lower.tail && !log.p) -dnorm((q-mean)/sd)*(D(x)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (lower.tail && log.p) dnorm((q-mean)/sd)*(D(x)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)/pnorm((q-mean)/sd)
  else -dnorm((q-mean)/sd)*(D(x)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)/pnorm((q-mean)/sd, lower.tail=FALSE))
  + stop("cannot take derivative wrt 'lower.tail' or 'log.p'")*(D(lower.tail) + D(log.p)),
  envir = sysDerivs)
newDeriv(dnorm(x, mean = 0, sd = 1, log = FALSE),
  (if (!log) dnorm((x-mean)/sd)/sd^2*((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1)) 
  else ((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1))/sd) 
  + stop("cannot take derivative wrt 'log'")*D(log),
  envir = sysDerivs)
newDeriv(asin(x), D(x)/sqrt(1+x^2), envir = sysDerivs)
newDeriv(acos(x), -D(x)/sqrt(1+x^2), envir = sysDerivs)
newDeriv(atan(x), D(x)/(1+x^2), envir = sysDerivs)
newDeriv(gamma(x), gamma(x)*digamma(x)*D(x), envir = sysDerivs)
newDeriv(lgamma(x), digamma(x)*D(x), envir = sysDerivs)
newDeriv(digamma(x), trigamma(x)*D(x), envir = sysDerivs)
newDeriv(trigamma(x), psigamma(x, 2L)*D(x), envir = sysDerivs)
newDeriv(psigamma(x, deriv = 0L), 
          psigamma(x, deriv + 1L)*D(x) 
        + stop("cannot take derivative wrt 'deriv'")*D(deriv), envir = sysDerivs)

newDeriv(x + y, D(x) + D(y), envir = sysDerivs)
newDeriv(x*y, x*D(y) + D(x)*y, envir = sysDerivs)
newDeriv(x/y, D(x)/y - x*D(y)/y^2, envir = sysDerivs)
newDeriv(x - y, D(x) - D(y), envir = sysDerivs)
newDeriv(x^y, y*x^(y-1)*D(x) + x^y*log(x)*D(y), envir = sysDerivs)
newDeriv((x), D(x), envir = sysDerivs)

# These are new

newDeriv(abs(x), sign(x)*D(x), envir = sysDerivs)
newDeriv(sign(x), 0, envir = sysDerivs)
