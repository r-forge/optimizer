# R-based replacement for deriv() function

sysDerivs <- new.env(parent = emptyenv())
sysSimplifications <- new.env(parent = emptyenv())

newDeriv <- function(expr, deriv, derivEnv = sysDerivs) {
    if (missing(expr))
    	return(ls(derivEnv))
    expr <- substitute(expr)
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    if (missing(deriv)) 
    	return(derivEnv[[fn]])
    deriv <- substitute(deriv)
    args <- expr[-1]
    argnames <- names(args)
    if (is.null(argnames))
    	argnames <- rep("", length(args))
    required <- which(argnames == "")
    for (i in required)
    	argnames[i] <- as.character(args[[i]])    	
    assign(fn, list(expr = expr, argnames = argnames, 
                    required = required, deriv = deriv), envir = derivEnv)
    invisible(derivEnv[[fn]])
}

newSimplification <- function(expr, test, simplification, do_eval = FALSE, simpEnv = sysSimplifications) {
    if (missing(expr))
    	return(ls(simpEnv))
    expr <- substitute(expr)
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    nargs <- length(expr) - 1L
    simps <- simpEnv[[fn]]
    if (missing(test)) {
        if (nargs <= length(simps))	
    	    return(simps[[nargs]])
    	else
    	    return(NULL)
    }
    test <- substitute(test)
    simplification <- substitute(simplification)
    
    args <- expr[-1]
    if (!is.null(names(args)))
    	stop("expr should not have named arguments")
    if (!all(sapply(args, is.name)))
    	stop("expr should have simple names as arguments")
    argnames <- sapply(args, as.character)
    if (any(duplicated(argnames)))
    	stop("expr names should be unique")
    	
    if (is.null(simps)) simps <- list()
    if (nargs <= length(simps)) 
    	simpn <- simps[[nargs]]
    else
    	simpn <- list()
    simpn <- c(simpn, list(list(expr = expr, argnames = argnames, test = test, 
                                simplification = simplification, do_eval = do_eval)))
    simps[[nargs]] <- simpn
    assign(fn, simps, envir = simpEnv)
}
    	
# This is a more general version of D()
Deriv <- function(expr, name, derivEnv = sysDerivs, do_substitute = TRUE, ...) {
    Recurse <- function(expr) {
    	if (is.call(expr)) {
    	    if (as.character(expr[[1]]) == "D")
    	    	expr <- Deriv(expr[[2]], name, derivEnv, do_substitute = FALSE, ...)
    	    else
    	    	for (i in seq_along(expr)[-1])
    	    	    expr[[i]] <- Recurse(expr[[i]])
    	}
    	expr
    }

    if (do_substitute)
    	expr <- substitute(expr)
    if (is.expression(expr))
    	return(as.expression(lapply(expr, Deriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, ...)))
    else if (is.numeric(expr) || is.logical(expr))
    	return(0)
    else if (is.call(expr)) {
    	fn <- as.character(expr[[1]])
	if (fn == "expression")
	    return(as.expression(lapply(as.list(expr)[-1], Deriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, ...)))
    	model <- derivEnv[[fn]]
    	if (is.null(model))
    	    stop("no derivative known for '", fn, "'")
 	if (missing(name)) {
 	    message("Pattern for")
 	    message(paste("  ", deparse(model$expr), collapse = "\n"))
 	    message("is")
 	    message(paste("  ", deparse(model$deriv), collapse = "\n"))
 	    return(invisible(NULL))
 	}
        args <- expr[-1]
        argnames <- names(args)
        if (is.null(argnames)) 
            argnames <- rep("", length(args))
        modelnames <- model$argnames
        argnum <- pmatch(argnames, modelnames)
        if (any(bad <- is.na(argnum) & argnames != ""))
            stop("Argument names not matched: ", paste(argnames[bad], collapse = ", "))
        unused <- setdiff(seq_along(modelnames), argnum)
        nonamecount <- sum(is.na(argnum))
        length(unused) <- nonamecount
        argnum[which(is.na(argnum))] <- unused
        default <- setdiff(seq_along(modelnames), argnum)
        if (length(bad <- setdiff(model$required, argnum)))
            stop("Missing required arguments: ", paste(modelnames[bad], collapse = ", "))
        
        # Now do the substitutions
        subst <- list()
        subst[argnum] <- as.list(args)
        subst[default] <- as.list(model$expr[-1])[default]
        names(subst) <- modelnames
        result <- do.call(substitute, list(model$deriv, subst))
        result <- Recurse(result)
        Simplify(result, ...)
    } else if (is.name(expr))
    	if (as.character(expr) == name)
    	    return(1)
    	else
    	    return(0)        
}

# This is a more general version of deriv(), since it allows user specified 
# derivatives and simplifications

fnDeriv <- function(expr, namevec, 
       hessian = FALSE, derivEnv = sysDerivs, ...) {
    if (!is.expression(expr)) expr <- as.expression(expr)
    if (length(expr) > 1)
	stop("Only single expressions allowed")
    exprs <- as.list(expr)
    n <- length(namevec)
    length(exprs) <- n + 1L
    for (i in seq_len(n))
	exprs[[i + 1]] <- Deriv(expr[[1]], namevec[i], derivEnv = derivEnv, do_substitute = FALSE, ...)
    names(exprs) <- c(".value", namevec)
    if (hessian) {
        m <- length(exprs)
	length(exprs) <- m + n*(n+1)/2
	for (i in seq_len(n))
	    for (j in seq_len(n-i+1) + i-1) {
		m <- m + 1
		exprs[[m]] <- Deriv(exprs[[i + 1]], namevec[j], derivEnv = derivEnv, 
		                    do_substitute = FALSE, ...)
	    }
    }
    exprs <- as.expression(exprs)
    subexprs <- findSubexprs(exprs)
    m <- length(subexprs)
    final <- subexprs[[m]]
    subexprs[[m]] <- substitute(.value <- expr, list(expr = final[[".value"]]))
    subexprs[[m+1]] <- substitute(.grad <- array(0, c(length(.value), 2L), list(NULL, namevec)),
				  list(namevec = namevec))
    if (hessian)
	subexprs[[m+2]] <- substitute(.hessian <- array(0, c(length(.value), 2L, 2L), 
					list(NULL, namevec, namevec)), 
					list(namevec = namevec))
    for (i in seq_len(n))
	subexprs[[m+1+i]] <- substitute(.grad[, name] <- expr,
			          list(name = namevec[i], expr = final[[namevec[i]]]))
    m <- length(subexprs)
    h <- 0
    if (hessian) {
	for (i in seq_len(n))
	    for (j in seq_len(n-i+1) + i-1) {
		h <- h + 1
		if (i == j)
		    subexprs[[m + h]] <- substitute(.hessian[, i, i] <- expr,
				  list(i = namevec[i], expr = final[[1 + n + h]]))
		else
		    subexprs[[m + h]] <- substitute(.hessian[, i, j] <- .hessian[, j, i] <- expr,
				  list(i = namevec[i], j = namevec[j], expr = final[[1 + n + h]]))
	    }	
	h <- h + 1
	subexprs[[m + h]] <- quote(attr(.value, "hessian") <- .hessian)
    }
    m <- length(subexprs)
    subexprs[[m+1]] <- quote(attr(.value, "gradient") <- .grad)
    subexprs[[m+2]] <- quote(.value)
    subexprs
}    

isFALSE <- function(x) identical(FALSE, x)
isZERO <- function(x) is.numeric(x) && length(x) == 1 && x == 0
isONE  <- function(x) is.numeric(x) && length(x) == 1 && x == 1
isMINUSONE <- function(x) is.numeric(x) && length(x) == 1 && x == -1

Simplify <- function(expr, simpEnv = sysSimplifications) {
    
    if (is.expression(expr))
    	return(as.expression(lapply(expr, Simplify, simpEnv)))    

    if (is.call(expr)) {
	for (i in seq_along(expr)[-1])
	    expr[[i]] <- Simplify(expr[[i]], simpEnv)
	fn <- as.character(expr[[1]])
	nargs <- length(expr) - 1
	while (!identical(simpEnv, emptyenv())) {
	    simps <- simpEnv[[fn]]
	    if (nargs > length(simps))
		return(expr)
	    simpn <- simps[[nargs]]
	    for (i in seq_along(simpn)) {  	
		argnames <- simpn[[i]]$argnames
		substitutions <- lapply(seq_along(argnames)+1L, function(i) expr[[i]])
		names(substitutions) <- argnames
		test <- simpn[[i]]$test
		if (with(substitutions, eval(test))) {
		    simplification <- do.call(substitute, list(simpn[[i]]$simplification, substitutions))
	   	    if (simpn[[i]]$do_eval)
	    		simplification <- eval(simplification)
	    	    return(simplification)
	        }
	     }
	     simpEnv <- parent.env(simpEnv)
	}
    }
    expr
}

findSubexprs <- function(expr, simplify = FALSE, tag = ".expr", ...) {
    digests <- new.env(parent = emptyenv())
    subexprs <- list()
    subcount <- 0
    
    record <- function(index) {
        if (simplify)
	    expr[[index]] <<- subexpr <- Simplify(expr[[index]], ...)
	else
	    subexpr <- expr[[index]]
	if (is.call(subexpr)) {
	    digest <- digest(subexpr)
	    for (i in seq_along(subexpr))
		record(c(index,i))
	    prev <- digests[[digest]]
	    if (is.null(prev)) 
		assign(digest, index, envir = digests)
	    else if (is.numeric(prev))  { # the index where we last saw this
	        subcount <<- subcount + 1
		name <- as.name(paste0(tag, subcount))
		assign(digest, name, envir = digests)
	    }
	}
    }
    
    edit <- function(index) {
	subexpr <- expr[[index]]
	if (is.call(subexpr)) {
	    digest <- digest(subexpr)	    
	    for (i in seq_along(subexpr))
		edit(c(index,i))
	    prev <- digests[[digest]]
	    if (is.name(prev)) {
		num <- as.integer(substring(as.character(prev), nchar(tag)+1L))
		subexprs[[num]] <<- call("<-", prev, expr[[index]])
		expr[[index]] <<- prev
	    } 
	}
    }

    
    for (i in seq_along(expr)) record(i)
    for (i in seq_along(expr)) edit(i)
    result <- quote({})
    result[seq_along(subexprs)+1] <- subexprs
    result[[length(result)+1]] <- expr
    result
}
    
# These are the derivatives supported by deriv()

newDeriv(log(x, base = exp(1)), 
         D(x)/(x*log(base)) - (1/base)*log(x, base)/log(base)*D(base))
newDeriv(exp(x), exp(x)*D(x))
newDeriv(sin(x), cos(x)*D(x))
newDeriv(cos(x), -sin(x)*D(x))
newDeriv(tan(x), 1/cos(x)^2*D(x))
newDeriv(sinh(x), cosh(x)*D(x))
newDeriv(cosh(x), sinh(x)*D(x))
newDeriv(sqrt(x), D(x)/2/sqrt(x))
newDeriv(pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),
  (if (lower.tail && !log.p) 
  	dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (!lower.tail && !log.p) 
  	-dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (lower.tail && log.p) 
  	(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, log = TRUE))
  else if (!lower.tail && log.p)	
  	-(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, lower.tail = FALSE, log = TRUE)))
  + stop("cannot take derivative wrt 'lower.tail' or 'log.p'")*(D(lower.tail) + D(log.p)))
newDeriv(dnorm(x, mean = 0, sd = 1, log = FALSE),
  (if (!log) 
  	dnorm((x-mean)/sd)/sd^2*((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1)) 
  else if (log)
        ((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1))/sd) 
  + stop("cannot take derivative wrt 'log'")*D(log))
newDeriv(asin(x), D(x)/sqrt(1+x^2))
newDeriv(acos(x), -D(x)/sqrt(1+x^2))
newDeriv(atan(x), D(x)/(1+x^2))
newDeriv(gamma(x), gamma(x)*digamma(x)*D(x))
newDeriv(lgamma(x), digamma(x)*D(x))
newDeriv(digamma(x), trigamma(x)*D(x))
newDeriv(trigamma(x), psigamma(x, 2L)*D(x))
newDeriv(psigamma(x, deriv = 0L), 
          psigamma(x, deriv + 1L)*D(x) 
        + stop("cannot take derivative wrt 'deriv'")*D(deriv))

newDeriv(x*y, x*D(y) + D(x)*y)
newDeriv(x/y, D(x)/y - x*D(y)/y^2)
newDeriv(x^y, y*x^(y-1)*D(x) + x^y*log(x)*D(y))
newDeriv((x), D(x))
# Need to be careful with unary + or -
newDeriv(`+`(x, y = .MissingVal), if (missing(y)) D(x) else D(x) + D(y))
newDeriv(`-`(x, y = .MissingVal), if (missing(y)) -D(x) else D(x) - D(y))

# These are new

newDeriv(abs(x), sign(x)*D(x))
newDeriv(sign(x), 0)

# Now, the simplifications

newSimplification(+a, TRUE, a)
newSimplification(-a, is.numeric(a), -a, do_eval = TRUE)
newSimplification(-a, is.call(a) && length(a) == 2 && as.character(a[[1]]) == "-", quote(a)[[2]], do_eval = TRUE)

newSimplification(exp(a), is.call(a) && length(a) == 2 && as.character(a[[1]]) == "log", quote(a)[[2]], do_eval = TRUE)
newSimplification(exp(a), is.numeric(a), exp(a), do_eval = TRUE)

newSimplification(log(a), is.call(a) && length(a) == 2 && as.character(a[[1]]) == "exp", quote(a)[[2]], do_eval = TRUE)
newSimplification(log(a), is.numeric(a), log(a), do_eval = TRUE)

newSimplification(!a, isTRUE(a), FALSE)
newSimplification(!a, isFALSE(a), TRUE)

newSimplification((a), TRUE, a)

newSimplification(a + b, isZERO(b), a)
newSimplification(a + b, isZERO(a), b)
newSimplification(a + b, identical(a, b), Simplify(quote(2*a)), do_eval = TRUE)
newSimplification(a + b, is.numeric(a) && is.numeric(b), a+b, do_eval = TRUE)

newSimplification(a - b, isZERO(b), a)
newSimplification(a - b, isZERO(a), Simplify(quote(-b)), do_eval = TRUE)
newSimplification(a - b, identical(a, b), 0)
newSimplification(a - b, is.numeric(a) && is.numeric(b), a - b, do_eval = TRUE)

newSimplification(a * b, isZERO(a), 0)
newSimplification(a * b, isZERO(b), 0)
newSimplification(a * b, isONE(a), b)
newSimplification(a * b, isONE(b), a)
newSimplification(a * b, isMINUSONE(a), Simplify(quote(-b)), do_eval = TRUE)
newSimplification(a * b, isMINUSONE(b), Simplify(quote(-a)), do_eval = TRUE)
newSimplification(a * b, is.numeric(a) && is.numeric(b), a * b, do_eval = TRUE)

newSimplification(a / b, isONE(b), a)
newSimplification(a / b, isMINUSONE(b), Simplify(quote(-a)), do_eval = TRUE)
newSimplification(a / b, isZERO(a), 0)
newSimplification(a / b, is.numeric(a) && is.numeric(b), a / b, do_eval = TRUE)

newSimplification(a ^ b, isONE(b), a)
newSimplification(a ^ b, is.numeric(a) && is.numeric(b), a ^ b, do_eval = TRUE)

newSimplification(log(a, base), is.call(a) && as.character(a[[1]]) == "exp", Simplify(call("/", quote(a)[[2]], quote(log(base)))), do_eval = TRUE)

newSimplification(a && b, isFALSE(a) || isFALSE(b), FALSE)
newSimplification(a && b, isTRUE(a), b)
newSimplification(a && b, isTRUE(b), a)

newSimplification(a || b, isTRUE(a) || isTRUE(b), TRUE)
newSimplification(a || b, isFALSE(a), b)
newSimplification(a || b, isFALSE(b), a)

newSimplification(if (cond) a, isTRUE(cond), a)
newSimplification(if (cond) a, isFALSE(cond), NULL)

newSimplification(if (cond) a else b, isTRUE(cond), a)
newSimplification(if (cond) a else b, isFALSE(cond), b)
newSimplification(if (cond) a else b, identical(a, b), a)

# This one is used to fix up the unary -
newSimplification(missing(a), identical(a, quote(.MissingVal)), TRUE)
newSimplification(missing(a), !identical(a, quote(.MissingVal)), FALSE)
