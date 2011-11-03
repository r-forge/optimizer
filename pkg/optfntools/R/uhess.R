############### uhess.R ####################
# ?? need tests of scaling to make sure we have everything right
uhess <- function(par, fnuser, ps=rep(1,length(par)), fs=1.0, maximize = FALSE, ...) {
    #         ihess<-ihess+1 ## possible counter
    # ?? attr(fnuser, counts)<-c(ifn, igr, ihess) ??? or attr(fnuser, ifn) etc.
    if (length(ps) == 1) ps<-rep(ps,length(par))
    npar <- length(par)
    th <- try(tryh <- fnuser$hess(par*ps, ...), silent = TRUE)
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
    if ((!is.null(maximize)) && maximize) 
        tryh <- -tryh # handle the maximization
    tryh<- diag(ps)%*%(tryh)%*%diag(ps)
    tryh<-tryh/fs
    attributes(tryh)<-tattr # restore attributes (?? not counting ihess??)
    tryh
}  # end uhess definition
############# end uhess ##########################
