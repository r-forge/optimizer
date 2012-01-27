############### ufn ####################
# function defined in order to deal with out of bounds functions/parameters
# ?? add exceeding function count inside and change attributes??  ##
ufn <- function(par, fnuser, ps=1.0, fs=1.0, maximize=FALSE, ...) {
#    if ((is.null(ps) || any(is.na(ps)) ) ) stop("No parameter scaling!")
    if (length(ps)>1 && any(is.na(ps)) ) stop("No parameter scaling!")
    else {
       ps<-rep(ps,length(par))
    }
    parps<-par*ps # Probably not necessary to pre-multiply
    testf <- try(tryf <- fnuser$fn(parps, ...), silent = TRUE)
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
    kkfn<-try(get("kfn",pos=sys.frame(fnuser$callpos)),silent=TRUE)
    if ((class(kkfn)!="try-error") && ! is.null(kkfn) && ! is.na(kkfn)) { 
       assign("kfn", kkfn+1, sys.frame(fnuser$callpos)) # pos=1 goes to .GlobalEnv
    } else  {
       warning("kfn NOT defined in calling program")
       assign("kfn", 1, sys.frame(fnuser$callpos))
    }
    if ((!is.null(maximize)) && maximize) tryf <- -tryf # handle the maximization
    tryf/fs # and scale to finish
}
# put probname into fnuser and change things that way
############## end ufn ###################

