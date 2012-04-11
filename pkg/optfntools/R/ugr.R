############### ugr.R ####################
ugr <- function(par, fnuser) {
    # Analytic gradient wrapper
    OPCON<-fnuser$OPCON
    npar <- length(par)
    if (is.null(fnuser$gr)){
       # Use numerical gradient
       deriv.approx<-attr(fnuser,"deriv.approx")
       if (is.null(deriv.approx)) deriv.approx=grfwd # default method is fwd diff
       if (is.null(fnuser$dots)) {
         tgr<-try(tryg<-deriv.approx(par*OPCON$PARSCALE, fnuser$fn), silent = TRUE)
       } else {
         tgr<-try(tryg<-deriv.approx(par*OPCON$PARSCALE, fnuser$fn, unlist(fnuser$dots)),
            silent = TRUE) 
       }
       #?? we are using the ORIGINAL function, scales parameters
    } else {
      if (is.null(fnuser$dots)) {
        tgr <- try(tryg <- fnuser$gr(par*OPCON$PARSCALE), silent = TRUE)
      } else {
        tgr <- try(tryg <- fnuser$gr(par*OPCON$PARSCALE, unlist(fnuser$dots)), 
               silent = TRUE)
      }
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
    OPCON$KGR<-1+OPCON$KGR
    if (OPCON$MAXIMIZE) tryg <- -tryg # Internal gradient
    tryg*OPCON$PARSCALE/OPCON$FNSCALE # handle the scaling
}
############# end ugr.R ##########################

