############### ugr.R ####################
ugr <- function(par, fnuser, ps=1.0, fs=1.0, maximize=FALSE, ...) {
    if (length(ps) == 1) ps<-rep(ps,length(par))
    # Analytic gradient wrapper
    # ?? add exceeding function count inside and change attributes?? #
    # ?? put numerical gradient inside this function ??#
    #   igr<-igr+1
    npar <- length(par)
    if (is.null(fnuser$gr)){
       # Use numerical gradient
       deriv.approx<-attr(fnuser,"deriv.approx")
       if (is.null(deriv.approx)) deriv.approx=grfwd # default method is fwd diff
       tgr<-try(tryg<-deriv.approx(par*ps, fnuser$fn, ...))
       #?? we are using the ORIGINAL function, scales parameters
    } else {
      tgr <- try(tryg <- fnuser$gr(par*ps, ...), silent = TRUE)
    }
    if ((class(tgr) == "try-error") || any(is.na(tryg)) || any(is.null(tryg)) || 
      any(is.infinite(tryg))) {
      tryg <- rep(.Machine$double.xmax, npar)
      attr(tryg, "inadmissible") <- TRUE
    }
    else {
      attr(tryg, "inadmissible") <- FALSE
    }
    if (any(is.null(tryg))) stop("NULL FUNCTION")
    kkgr<-try(get("kgr",pos=sys.frame(fnuser$callpos)),silent=TRUE)
    if ((class(kkgr)!="try-error") && ! is.null(kkgr) && ! is.na(kkgr)) { 
       assign("kgr", kkgr+1, pos=sys.frame(fnuser$callpos))
    } else {
       warning("kfn undefined in calling program -- set to 1")
       assign("kfn", 1, sys.frame(fnuser$callpos))
    }
    if ((!is.null(maximize)) && maximize) 
      tryg <- -tryg # Internal gradient
    tryg*ps/fs # handle the scaling
}
############# end ugr.R ##########################

