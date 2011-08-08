############### ufn ####################
# function defined in order to deal with out of bounds functions/parameters
# ?? add exceeding function count inside and change attributes??  ##
ufn <- function(par, fnuser, ps=rep(1.0, length(par)), fs=1.0, maximize = FALSE, ...) {
    if (length(ps) == 1) ps<-rep(ps,length(par))
    testf <- try(tryf <- fnuser$fn(par*ps, ...), silent = TRUE)
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
