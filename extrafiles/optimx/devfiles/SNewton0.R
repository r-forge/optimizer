SNewton <- function(x0, fn, gr, hess, variant=list(ls<-1, solve=1, marq=1), 
                    control = list(trace = 2, maxit = 1000), ...) {
    ## Safeguarded Newton minimizer
    ##
    ##Input
    ##       - x0 is the initial value
    ##       - fn is the function we wish to minimize
    ##       - gr is its gradient function
    #            Probably do NOT want to allow numericals
    ##       - hess is its hessian 
    #            Need to allow to be NULL and use numericals
    #            via Jacobian
    ##       - variant is list specifying line search, solver and marquardt strategy
    #  Possibly spec marquardt via lambda and phi parameters -- 0 means no marquardt
    ##       - ... is data used in the function fn
    ##Output (list) -- need to match optim() output!! ???
    ##       - xs is the value at the minimum
    ##       - fv is the fn evaluated at xs
    ##       - grd is the gradient
    ##       - Hess is the Hessian
    ##       - niter is the number of interations needed.
    
    ## counters
    nfn<-0 # function evaluations
    ngr<-0 # gradient evaluations
    nhess<-0 # hessian evaluations
    niter<-1 # Newton iterations

    ## Should ensure strategy is included in output in words ??


    # -------- Line Search approach 1 ----------  !! CLEAN UP !!
    lnsrch <- function(fn, lgr=NULL, xc, d, ...) { # line search approach 1
        ## put in lgr to allow for gradient in line search
        ## Uses Brent's method to find the best stepsize gamma in [0.1,1]
        flsch <- function(gm, fn, xc, d, ...) {
            fval <- fn(xc + gm * d, ...)
            fval
        }
        lout <- optimize(flsch, interval = c(0.1, 1), fn = fn, 
            xc = xc, d = d, ...)$min
        lout
    }
    # ---- end line search approach 1 ----
    # -------- Line Search approach 2 ----------  !! CLEAN UP !!
    # Use slope + fval0 + fvalstep and interpolate as in Rcgmin (??)
    # ---- end line search approach 2 ----
    # -------- Line Search approach 3 ----------  !! CLEAN UP !!
    # Use backtrack as in Rvmmin (??)
    # ---- end line search approach 3 ----
    
    ## CONTROLS:
    # ?? may want to check maxit etc. and other values for sense
    eps0 <- .Machine$double.eps
    eps <- (eps0)^(1/3)
    lambda<-0.0001
    decr<-0.4 # Factor to decrease Levenberg Marquardt lambda
    incr<-10 # Factor to increase Levenberg Marquardt lambda
    ## ---- end controls ----
    if (control$trace > 0) {
	cat("Safeguarded Newton method for function minimization\n")
        cat("J C Nash 2010")
        # ??variant
    }
    # Compute base function
    f0<-fn(x0, ...)
    nfn<-nfn+1
    while (niter <= control$maxit) {
        grd <- gr(x0, ...)
        
        H <- hess(x0, ...)
        Haug <- H + lambda*(diag(H) + lambda) * eps
        stp <- solve(H, -grd)
        ## Do line search
        gvl <- lnsrch(fn, x0, stp, ...)
        xn <- c(x0 + gvl * stp)
        if (max(abs(grd)) < eps) 
            break
        if (niter == control$maxit) {
            cat("NewtonR: Failed to converge!\n")
            return(0)
        }
        fnn<-fn(xn,...)
        cat(niter,"  f=",fnn,"\n")
        x0 <- xn
        niter<-1+niter
    }
    # ?? work out "convergence" message
    out <- NULL
    out$xs <- xn
    out$fv <- fn(xn, ...)
    out$grd <- grd
    out$Hess <- H
    out$niter <- niter
    out
}




