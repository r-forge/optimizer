############### ugr.R ####################
ugr <- function(par, fnuser, maximize=FALSE, ...) {
    # Analytic gradient wrapper
    # ?? add exceeding function count inside and change attributes?? #
    #   igr<-igr+1
    npar <- length(par)
    tgr <- try(tryg <- fnuser$gr(par, ...), silent = TRUE)
    if ((class(tgr) == "try-error") || any(is.na(tryg)) || any(is.null(tryg)) || 
        any(is.infinite(tryg))) {
        tryg <- rep(.Machine$double.xmax, npar)
        attr(tryg, "inadmissible") <- TRUE
    }
    else {
        attr(tryg, "inadmissible") <- FALSE
    }
    if (any(is.null(tryg))) 
        stop("NULL FUNCTION")
    #         attr(tryg,'igr')<-igr
    if ((!is.null(maximize)) && maximize) 
        tryg <- -tryg
    # handle the maximization
    tryg
}
############# end ugr ##########################
