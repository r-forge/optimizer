############### ufn ####################
# function defined in order to deal with out of bounds
#   functions/parameters
# ?? add exceeding function count inside and change attributes??  ##
#   nfun<-nfun+1
ufn <- function(par, fnuser, maximize = FALSE, ...) {
    testf <- try(tryf <- fnuser$fn(par, ...), silent = TRUE)
    # try to Compute the function. Should we quote it?
    if ((class(testf) == "try-error") || is.na(tryf) || is.null(tryf) || 
        is.infinite(tryf)) {
        tryf <- .Machine$double.xmax
        attr(tryf, "inadmissible") <- TRUE
    }
    else {
        attr(tryf, "inadmissible") <- FALSE
    }
    if (is.null(tryf)) 
        stop("NULL FUNCTION")
    #    attr(tryf, 'nfun')<-nfun
    if ((!is.null(maximize)) && maximize) 
        tryf<- -tryf # handle the maximization
    tryf
}
############## end ufn ###################

