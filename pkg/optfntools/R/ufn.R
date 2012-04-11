############### ufn ####################
# function defined in order to deal with out of bounds functions/parameters
ufn <- function(par, fnuser) {
    OPCON<-fnuser$OPCON
    parps<-par*OPCON$PARSCALE 
    if (is.null(fnuser$dots)) {
       testf <- try(tryf <- fnuser$fn(parps), silent = TRUE)
    } else {
       testf <- try(tryf <- fnuser$fn(parps, unlist(fnuser$dots)), silent = TRUE)
    }
    if ((class(testf) == "try-error") || is.na(tryf) || is.null(tryf) || 
        is.infinite(tryf)) {
        tryf <- .Machine$double.xmax
        attr(tryf, "inadmissible") <- TRUE
    }
    else {
        attr(tryf, "inadmissible") <- FALSE
    }
    if (is.null(tryf)) stop("NULL FUNCTION")
    OPCON$KFN<-1+OPCON$KFN
    if (OPCON$MAXIMIZE) tryf <- -tryf # handle the maximization
    tryf/OPCON$FNSCALE # and scale to finish
}
# put probname into fnuser and change things that way
############## end ufn ###################


