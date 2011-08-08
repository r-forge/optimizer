############### uhess.R ####################
# ?? need tests of scaling to make sure we have everything right
uhess <- function(par, fnuser, ps=rep(1,length(par)), fs=1.0, maximize = FALSE, ...) {
    #         ihess<-ihess+1 ## possible counter
    if (length(ps) == 1) ps<-rep(ps,length(par))
    npar <- length(par)
    th <- try(tryh <- fnuser$hess(par*ps, ...), silent = TRUE)
    if ((class(th) == "try-error") || any(is.na(tryh)) || any(is.null(tryh)) || 
        any(is.infinite(tryh))) {
        tryh <- matrix(.Machine$double.xmax, nrow = npar, ncol = npar)
        attr(tryh, "inadmissible") <- TRUE
    }
    else {
        attr(tryh, "inadmissible") <- FALSE
    }
    if (any(is.null(tryh))) 
        stop("NULL FUNCTION")
    #         attr(tryh,'ihess')<-ihess
    if ((!is.null(maximize)) && maximize) 
        tryh <- -tryh # handle the maximization
     tryh<- diag(ps)%*%(tryh)%*%diag(ps)
     tryh/fs
}  # end uhess definition
############# end uhess ##########################
