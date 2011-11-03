# fakeopt is below -- here are the useroptfn functions

debug<-FALSE


############### ufn ####################
# function defined in order to deal with out of bounds
#   functions/parameters
# ?? add exceeding function count inside and change attributes??  ##
#   nfun<-nfun+1
#
#      Gradient: usergr is fn of parmuser 
#                mgr is usergr(parmuser) * ps / fs
#            =   d (mfn) / d(mparm) = (d (userfn/fs) / d parmuser) * (d parmuser / d mparm)
#            =   usergr * (1/fs) * ps
#
ufn <- function(par, fnuser, ps=rep(1.0, length(par)), fs=1.0, maximize = FALSE, ...) {
    if (length(ps) == 1) ps<-rep(ps,length(par))
    if (debug) {
    cat("In ufn, fs=", fs, "\n")
    cat("   ps:")
    print(ps)
    cat("maximize=",maximize,"\n")
    cat("raw parameters:")
    print(par)
    cat(" scaled parameters:")
    print(par*ps)
    fargsx<-list(...)
    cat("fargsx:")
    print(fargsx)
    cat("fnuser$fn")
    print(fnuser$fn)
    cat("about to compute tryf\n")
    tryf <- fnuser$fn(par*ps, ...)
    cat("First try at function =",tryf,"\n")
    cat("Params to user fn:")
    print(par*ps)
    }
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
    if (debug) {
    cat("maximize:")
    print(maximize)
    }
    if ((!is.null(maximize)) && maximize) tryf <- -tryf # handle the maximization
    if (debug) {
    cat("unscaled function is ",tryf,"  finishing ufn\n")
    }
    tryf/fs # and scale to finish
}
############## end ufn ###################

############### ugr.R ####################
ugr <- function(par, fnuser, ps=rep(1,length(par)), fs=1.0, maximize = FALSE, ...) {
    if (length(ps) == 1) ps<-rep(ps,length(par))
    # Analytic gradient wrapper
    # ?? add exceeding function count inside and change attributes?? #
    #   igr<-igr+1
    npar <- length(par)
    tgr <- try(tryg <- fnuser$gr(par*ps, ...), silent = TRUE)
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
        tryg <- -tryg # Internal gradient
    tryg*ps/fs # handle the scaling
}
############# end ugr ##########################

############### uhess.R ####################
# ?? need tests of scaling to make sure we have everything right
uhess <- function(par, fnuser, ps=rep(1,length(par)), fs=1.0, maximize = FALSE, ...) {
    #         ihess<-ihess+1 ## possible counter
    if (length(ps) == 1) ps<-rep(ps,length(par))
    npar <- length(par)
    print(ps)
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
     print(tryh)
#    tryh*ps*ps/fs # handle the scaling ?? check
     tryh<- diag(ps)%*%(tryh)%*%diag(ps)
     tryh/fs
}  # end uhess definition
############# end uhess ##########################




fakeopt <- function(par, fn, gr = NULL, lower = NULL, 
    upper = NULL, bdmsk = NULL, control = list(), ...) {
    # Input:
    #  par = a vector containing the starting point
    # fn = objective function (assumed to be sufficeintly
    #   differentiable)
    #  gr = gradient of objective function
    #  lower = vector of lower bounds on parameters
    #  upper = vector of upper bounds on parameters
    # Note: free parameters outside bounds will be adjusted to
    #   bounds.
    # bdmsk = control vector for bounds and masks. Parameters
    #   for which bdmsk are 1
    # are unconstrained or 'free', those with bdmsk 0 are
    #   masked i.e., fixed.
    #  control = list of control parameters
    #  ... = dotargs -- extras
    ##
    # control defaults -- taken from spg -- need to adjust for
    #   Rcgmin and Rvmmin
    # NOT yet in control set
    dblmax<-.Machine$double.xmax
    keepinputpar <- TRUE  # Do not let bmchk change parameters to nearest bound
    #  ?? put keepinputpar into controls??
    ctrl <- list(maxit = 500, maxfeval = 3000, maximize = FALSE, 
        trace = 0, eps = 1e-07, usenumDeriv = FALSE, dowarn = TRUE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control  #
    #  M           <- ctrl$M      # not needed
    maxit <- ctrl$maxit  #
    maxfeval <- ctrl$maxfeval  #
    maximize <- ctrl$maximize  # TRUE to maximize the function
    trace <- ctrl$trace  #
    eps <- ctrl$eps  #
    dowarn <- ctrl$dowarn  #
    grNULL <- is.null(gr)  # if gr function is not provided, we want to use numDeriv
    fargs <- list(...)  # the ... arguments that are extra function / gradient data
    #################################################################
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) 
        cat("fakeopt -- to test user functions\n")
    bvec <- par  # copy the parameter vector
    n <- length(bvec)  # number of elements in par vector
    if (trace > 0) {
        cat("Problem of size n=", n, "  Extra arguments:\n")
        print(fargs)
    }
    ifn <- 1  # count function evaluations
   btest <- bmchk(bvec, lower = lower, upper = upper, bdmsk = bdmsk, 
        trace = trace)
   cat("Results of bmchk: ")
   print(btest)
   tmp<-readline("continue")
   bvec<-bvec/fargs$ps
   if(trace>0) {
      cat("scaled bvec/ps:")
      print(bvec)
   }
   if (!btest$admissible) 
        stop("Inadmissible bounds: one or more lower>upper")
    if (btest$parchanged) {
        if (keepinputpar) 
            stop("Parameter out of bounds")
        else warning("Parameter out of bounds has been moved to nearest bound")
    }
    nolower <- btest$nolower
    noupper <- btest$noupper
    bounds <- btest$bounds
    bdmsk <- btest$bdmsk  # change bdmsk to values set in bmchk
    if (trace > 3) {
        cat("Adjusted bdmsk vector:")
        print(bdmsk)
    }
    lower <- btest$lower
    upper <- btest$upper
    ############## end bounds check #############
    f <- try(do.call("fn", append(list(bvec), fargs)), silent = TRUE)  # Compute the function.
    if ((class(f) == "try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
        cat("FAILED!\n")
        msg <- "Initial point gives inadmissible function value"
        conv <- 20
        if (trace > 0) 
            cat(msg, "\n")
        # change NA to dblmax 110524
        ans <- list(bvec, dblmax, c(ifn, 0), 0, conv, msg, bdmsk)  #
        names(ans) <- c("par", "value", "counts", "convergence", 
            "message", "bdmsk")
        return(ans)
    }
    if (trace > 0) 
        cat("Initial fn=", f, "\n")
    if (trace > 2) 
        print(bvec)
    keepgoing <- TRUE  # to ensure loop continues until we are finished
    ig <- 1  # count gradient evaluations
    ilast <- ig  # last time we used gradient as search direction
    fmin <- f  # needed for numerical gradients
# ?? need to ensure we pass through the fnuser
    g <- ugr(bvec, ...)  # Do we need to use try() ?? Possible not
    if (trace > 2) {
        cat("g:")
        print(g)
    }
# --- we don't actually minimize here
    if (trace > 0) 
        cat("Seem to be done fakeopt\n")
    if (maximize) 
        fmin <- (-1) * fmin
    conv<-2
    msg<-"this is the end of fakeopt"
    ans <- list(par, fmin, c(ifn, ig), conv, msg, bdmsk)
    ## ?? need to fix up message
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message", "bdmsk")
    #return(ans)
    ans
}  ## end of fakeopt
# --- end fakeopt ---

