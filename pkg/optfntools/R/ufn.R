############### ufn ####################
# function defined in order to deal with out of bounds functions/parameters
# ?? add exceeding function count inside and change attributes??  ##
ufn <- function(par, fnuser, ps=1.0, fs=1.0, maximize=FALSE, ...) {
    if (any(is.na(ps))) stop("No parameter scaling!")
    if (length(ps) == 1) {
       ps<-rep(ps,length(par))
    }
    parps<-par*ps # Probably not necessary to pre-multiply
    testf <- try(tryf <- fnuser$fn(parps, ...))
    # try to Compute the function. Should we quote it?
    if ((class(testf) == "try-error") || is.na(tryf) || is.null(tryf) || 
        is.infinite(tryf)) {
        tryf <- .Machine$double.xmax
        attr(tryf, "inadmissible") <- TRUE
    }
    else {
        attr(tryf, "inadmissible") <- FALSE
    }
    if (is.null(tryf)) stop("NULL FUNCTION")
    #    attr(tryf, 'nfun')<-nfun
    if ((!is.null(maximize)) && maximize) tryf <- -tryf # handle the maximization
    tryf/fs # and scale to finish
}
############## end ufn ###################
