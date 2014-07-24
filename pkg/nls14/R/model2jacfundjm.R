model2resfundjm <- function(modelformula, pvec, data = NULL, testresult = TRUE) {
    stopifnot(inherits(modelformula, "formula"))

    if (length(modelformula) == 2) {
        residexpr <- modelformula[[2]]
    } else if (length(modelformula) == 3) {
        residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    } else stop("Unrecognized formula")
    
    if (is.null(names(pvec)))
	names(pvec) <- paste0("p", seq_along(pvec))
    
    if (is.null(data))
	data <- environment(modelformula)
    else if (is.list(data))
	data <- list2env(data, parent = environment(modelformula))
    else if (!is.environment(data))
	stop("'data' must be a dataframe, list, or environment")
    
    mfun <- function(prm) {
        if (is.null(names(prm))) 
	    names(prm) <- names(pvec)
	localdata <- list2env(as.list(prm), parent = data)
	eval(residexpr, envir = localdata)
    }
    
    if (testresult) {
	resids <- mfun(pvec)
	if (any(!is.finite(resids))) 
	    stop("residuals contain ", unique(resids[!is.finite(resids)]))
	rm(resids)
    }
    
    mfun
}

model2jacfundjm <- function(modelformula, pvec, data = NULL, testresult = TRUE) {

    stopifnot(inherits(modelformula, "formula"))

    if (length(modelformula) == 2) {
        residexpr <- modelformula[[2]]
    } else if (length(modelformula) == 3) {
        residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    } else stop("Unrecognized formula")
    
    if (is.null(names(pvec)))
	names(pvec) <- paste0("p", seq_along(pvec))
    
    derivexpr <- deriv(residexpr, names(pvec))
    
    if (is.null(data))
	data <- environment(modelformula)
    else if (is.list(data))
	data <- list2env(data, parent = environment(modelformula))
    else if (!is.environment(data))
	stop("'data' must be a dataframe, list, or environment")
    
    jfun <- function(prm) {
        if (is.null(names(prm))) 
	    names(prm) <- names(pvec)
	localdata <- list2env(as.list(prm), parent = data)
	attr(eval(derivexpr, envir = localdata), "gradient")
    }
    
    if (testresult) {
	derivs <- jfun(pvec)
	if (any(!is.finite(derivs))) 
	    stop("derivatives contain ", unique(derivs[!is.finite(derivs)]))
	rm(derivs)
    }
    
    jfun
}
