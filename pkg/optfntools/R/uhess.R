############### uhess.R ####################
# ?? need tests of scaling to make sure we have everything right
uhess <- function(par, fnuser) {
    OPCON<-fnuser$OPCON
    npar <- length(par)
    if (is.null(fnuser$dots)) {
       th <- try(tryh <- fnuser$hess(par*OPCON$PARSCALE), silent = TRUE)
    } else {
       th <- try(tryh <- fnuser$hess(par*OPCON$PARSCALE, 
             unlist(fnuser$dots)), silent = TRUE)
    }
    if ((class(th) == "try-error") || any(is.na(tryh)) || any(is.null(tryh)) || 
        any(is.infinite(tryh))) {
        tryh <- matrix(.Machine$double.xmax, nrow = npar, ncol = npar)
        attr(tryh, "inadmissible") <- TRUE
        cat("INADMISSIBLE\n")
    } else {
        attr(tryh, "inadmissible") <- FALSE # assume OK until proven otherwise
    }
    tattr<-attributes(tryh) # save all attrubutes
    if (any(is.null(tryh))) stop("NULL FUNCTION")
    OPCON$KHESS<-1+OPCON$KHESS
    if (OPCON$MAXIMIZE) tryh <- -tryh # handle the maximization
    tryh<- diag(OPCON$PARSCALE)%*%(tryh)%*%diag(OPCON$PARSCALE) 
    # attributes NOT inherited in operation above
    tryh<-tryh/OPCON$FNSCALE
    attributes(tryh)<-tattr # restore attributes (?? not counting ihess??)
    tryh
}  # end uhess definition
############# end uhess ##########################
