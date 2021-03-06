optimr <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {


  npar <- length(par)
  defctrl <- ctrldefault(npar) # could leave this out in most cases
  if (is.null(method)) method <- defctrl$defmethod

  outmethod <- checksolver(method, defctrl$allmeth, defctrl$allpkg) # there will only be one! 
  if (is.null(outmethod)) {
		if (control$trace > 0) cat("Solver ",method," missing\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 8888 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- paste("Missing method ",method)
                ans$hessian <- NULL
                return(ans) # can't proceed without solver
  }



# Check if bounded
  bdmsk <- bmchk(par, lower=lower, upper=upper)
  control$have.bounds <- bdmsk$bounds # and set a control value

  orig.method <- method
  orig.gr <- gr
  orig.fn <- fn
  savehess <- hessian


  if (is.null(control$trace)) control$trace <- defctrl$trace

  if (is.null(control$parscale)) { 
        pscale <- rep(1,npar)
        if(control$trace > 0) { cat("Unit parameter scaling\n") }
  } else { 
        pscale <- control$parscale 
        if(control$trace > 0) {
          cat("Parameter scaling:")
          print(pscale)
        }
  }
  spar <- par/pscale # scaled parameters
  slower <- -Inf
  supper <- Inf # to ensure defined
  if (control$have.bounds) {
    slower <- lower/pscale
    supper <- upper/pscale
  }
  fnscale <- 1 # default to ensure defined
  if (is.null(control$fnscale)) {
     if (! is.null(control$maximize) && control$maximize ) {fnscale <- -1}
  else if (! is.null(control$maximize)) {
          if ( (control$fnscale < 0) && control$maximize) {fnscale <- -1} # this is OK
          else stop("control$fnscale and control$maximize conflict")
       } # end ifelse
  } # end else
  control$fnscale <- fnscale # to ensure set again

# 160615 -- decided to abandon nloptr in optimz

  efn <- function(spar, ...) {
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn(par, ...) * fnscale
  }

  appgr<-FALSE # so far assuming analytic gradient
#  DO NOT PROVIDE A DEFAULT -- LET METHOD DO THIS
#  if (is.null(gr)) gr <- defctrl$defgrapprox
#  if (is.null(gr)) cat("gr is NULL\n")
  if (is.character(gr)) {
     appgr <- TRUE # to inform us that we are using approximation
     egr <- function(spar, ...){
        if (control$trace > 1) {
           cat("fnscale =",fnscale,"  pscale=")
           print(pscale)
           cat("gr:")
           print(gr)
           cat("par:")
           print(par)
        }
        par <- spar*pscale
        result <- do.call(gr, list(par, userfn=fn, ...)) * fnscale
     }
  } else { 
    if (is.null(gr)) {egr <- NULL}
    else {
       egr <- function(spar, ...) {
         par <- spar*pscale
         result <- gr(par, ...) * pscale * fnscale
       }
    }
  } # end egr definition

  if (is.null(hess)) { ehess <- NULL}
  else { ehess <- function(spar, ...) {
                      par <- spar*pscale
                      result <- hess(par, ...) * pscale * pscale * fnscale
                      result
## ?? NEED TO CHECK THIS
                  }
  }

  if (appgr && (control$trace>0)) cat("Using numerical approximation '",gr,"' to gradient in optimru()\n")

  nlmfn <- function(spar, ...){
     f <- efn(spar, ...)
     if (is.null(egr)) {g <- NULL} else {g <- egr(spar, ...)}
     attr(f,"gradient") <- g
     if (is.null(ehess)) { h <- NULL } else {h <- ehess(spar, ...)}
     attr(f,"hessian") <- h
     f
  }
# ?? do we want ehess ?    Not at 150714

## Masks 
   maskmeth <- defctrl$maskmeth
   msk <- bdmsk$bdmsk # Only need the masks bit from here on
   if (any(msk == 0) ) {
      if ( !(method %in% maskmeth) ) {
         stopmsg <- paste("Method ",method," cannot handle masked (fixed) parameters")
         stop(stopmsg)
      }
      if (control$trace > 0) cat("Masks present\n")
   }

# replacement for optim to minimize using a single method

# time is in opm(), but not here
# The structure has   par, value, counts, convergence, message, hessian

# Run a single method

# expand bounds
  if (length(lower) == 1) lower<-rep(lower,npar)
  if (length(upper) == 1) upper<-rep(upper,npar)

  mcontrol <- list() # define the control list

# Methods from optim()
  if (method == "Nelder-Mead" || 
      method == "BFGS" || 
      method == "L-BFGS-B" || 
      method == "CG" || 
      method == "SANN") {
      # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
      mcontrol$maxit <- defctrl$maxit # 160922 
      mcontrol$trace <- control$trace
##	mcontrol$parscale <- control$parscale # Use internal scaling
      mcontrol$parscale <- NULL # using user fn 

# Note: hessian always FALSE

#        cat("Before optim() call - control$have.bounds =",control$have.bounds,"\n")
      if (control$have.bounds) {
        if (method != "L-BFGS-B") {
            errmsg <- "optim() can only handle bounds with L-BFGS-B\n"
            if (control$trace > 0) cat(errmsg,"\n")
            ans <- list()
            class(ans)[1] <- "try-error"
            warning("optimr: optim() with bounds ONLY uses L-BFGS-B")
        } else {
#            ans <- try(optim(par=par, fn=orig.fn, gr=orig.gr, 
#                      lower=lower, upper=upper, method="L-BFGS-B", hessian=FALSE, 
#                       control=mcontrol, ...))
            ans <- try(optim(par=par, fn=efn, gr=egr, 
                      lower=lower, upper=upper, method="L-BFGS-B", hessian=FALSE, 
                       control=mcontrol, ...))
          }
        } else {
#          cat("calling optim() with no bounds\n")
#          ans <- try(optim(par=par, fn=orig.fn, gr=orig.gr, 
#                method=method, hessian=FALSE, control=mcontrol, ...))
          ans <- try(optim(par=par, fn=efn, gr=egr, 
                method=method, hessian=FALSE, control=mcontrol, ...))
#          print(ans)

        }
        if (class(ans)[1] == "try-error") { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
                errmsg <- "optim method failure\n"
                if (method != "L-BFGS-B") errmsg <- paste("optim() with bounds ONLY uses L-BFGS-B: ", errmsg)
		if (control$trace>0) cat(errmsg)
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- errmsg
        } # otherwise ans is OK and we return it
        ## return(ans) # to ensure we return
      }   # end if using optim() methods
## --------------------------------------------
      else if (method == "nlminb") {
        # Here we use portLib routine nlminb rather than optim as our minimizer
        mcontrol$iter.max<-mcontrol$maxit # different name for iteration limit in this routine
        mcontrol$maxit<-NULL # and we null it out
        mcontrol$abs.tol <- 0 # To fix issues when minimum is less than 0. 20100711
        mcontrol$eval.max <- control$maxfeval
	if ( is.null(control$trace) || is.na(control$trace) || control$trace == 0) { 
		mcontrol$trace = 0
	} else { 
		mcontrol$trace = 1 # this is EVERY iteration. nlminb trace is freq of reporting.
	}
        ans <- try(nlminb(start=spar, objective=efn, gradient=egr, lower=slower, 
		upper=supper, control=mcontrol,  ...))
        if (class(ans)[1] != "try-error") {
		# Translate output to common format and names
        	ans$value<-ans$objective
                ans$par <- ans$par*pscale
	        ans$objective<-NULL
	        ans$counts[1] <- ans$evaluations[1]
        	ans$counts[2] <- ans$evaluations[2]
		ans$evaluations<-NULL # cleanup
	        ans$iterations<-NULL
                ans$hessian <- NULL
	} else { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		if (control$trace>0) cat("nlminb failure\n")
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using nlminb
## --------------------------------------------
      else if (method == "nlm") { # Use stats package nlm routine
#        if (is.null(gr)) { stop("optimr -- nlm -- we do not allow gr = NULL") }
	if (! is.null(control$maxit) ) {iterlim <- control$maxit }
        else { iterlim <- 100 }
	print.level <- 0 
        errmsg <- NULL
        if (control$have.bounds) {
              if(control$trace > 0) cat("nlm cannot handle bounds\n")
              errmsg <- "nlm cannot handle bounds\n"
            ##  stop("nlm tried with bounds")
            ans <- list()
            class(ans)[1] <- "try-error"
        } else {
          if (! is.null(control$trace) && (control$trace > 0) ) {print.level <- 2 } 
          ans <- try(nlm(f=nlmfn, p=spar, iterlim=iterlim, print.level=print.level, ...))
        }
        if (class(ans)[1] != "try-error") {
		if (ans$code == 1 || ans$code == 2 || ans$code == 3) ans$convergence <- 0
		if (ans$code == 4) ans$convergence <- 1
                if (ans$code == 5) ans$convergence <- 5
        	# Translate output to common format
		ans$value <- ans$minimum
		ans$minimum <- NULL
                ans$par <- ans$estimate*pscale
		ans$estimate <- NULL
        	ans$counts[2] <- ans$iterations
                ans$counts[1] <- NA
        	ans$iterations <- NULL
                ans$hessian <- NULL
                ans$gradient <- NULL # We lose information here
                ans$message <- paste("Convergence indicator (code) = ",ans$code)
                ans$code <- NULL
	} else {
		if (control$trace > 0) cat("nlm failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        print.level <- NULL # clean up
        ## return(ans)
      } # end if using nlm
## --------------------------------------------
      else if (method == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
        mcontrol$trace <- control$trace
        mcontrol$maxit <- control$maxit # 151217 JN
        if (! is.null(egr)) {
  	  if (control$have.bounds) { # 151220 -- this was not defined
            # 170919 -- explicit reference to package
   	    ans <- try(Rcgmin::Rcgminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
	  } else {
   	     ans <- try(Rcgmin::Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
	  }
        }
        if (!is.null(egr) && (class(ans)[1] != "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (control$trace>0) {
                    cat("Rcgmin failed for current problem \n")
                    if(is.null(egr)) cat("Note: Rcgmin needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                if(is.null(egr)) {
                   ans$message <- "Must specify gradient function for Rcgmin"       
                   ans$convergence <- 9998 # for no gradient where needed
                }
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (method == "Rvmmin") { # Use Rvmmin routine (ignoring masks??)
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
	mcontrol$trace <- control$trace # 140902 Note no check on validity of values
	if (! is.null(egr)) {
          if(control$have.bounds) { # 170919 make package explicit
   	     ans <- try(Rvmmin::Rvmminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
	  } else {
             ans <- try(Rvmmin::Rvmminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
	  }
        }
        if (! is.null(egr) && (class(ans)[1] != "try-error")) {
            ans$par <- ans$par*pscale
            ans$bdmsk <- NULL
        } else {
            if (control$trace>0) {
                cat("Rvmmin failed for current problem \n")
                if(is.null(egr)) cat("Note: Rcgmin needs gradient function specified\n")
            }
	    ans<-list() # ans not yet defined, so set as list
            ans$convergence <- 9999 # failed in run
	    ans$value <- defctrl$badval
	    ans$par<-rep(NA,npar)
	    ans$counts[1] <- NA # save function and gradient count information
	    ans$counts[2] <- NA # save function and gradient count information
	    ans$message <- NULL        
            if(is.null(egr)) {
               ans$message <- "Must specify gradient function for Rvmmin"
               ans$convergence <- 9998 # for no gradient where needed
            }
            ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rvmmin
## --------------------------------------------
      else if (method == "hjn") {# Use JN Hooke and Jeeves
        if (control$trace > 0) { # this function is in optimr, so does not need explicit package
           cat("hjn:control$have.bounds =",control$have.bounds,"\n")
           cat("optimr - hjn - msk:")
           print(msk)
        }
        ans <- try(hjn(spar, efn, lower=slower, upper=supper, bdmsk=msk, 
                        control=control, ...))
        if (class(ans)[1] != "try-error") {
            ## Need to check these carefully??
            ans$par <- ans$par*pscale
            ans$value <- ans$value*fnscale
            ans$message <- NA # Should add a msg ??
         } else {
            if (control$trace > 0) cat("hjn failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- defctrl$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         ## return(ans)
      }  ## end if using hjn
## --------------------------------------------
      else if (method == "spg") { # Use BB package routine spg as minimizer
        mcontrol$maximize <- NULL # Use external maximization approach
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
        if (control$trace > 0) { 
            mcontrol$trace <- TRUE
            if (control$trace > 1) mcontrol$triter <- 1 # default is 10
        } else { mcontrol$trace <- FALSE }
        ans <- try(BB::spg(par=spar, fn=efn, gr=egr, lower=slower, upper=supper,  
		control=mcontrol, ...))
        if (class(ans)[1] != "try-error") { 
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval<-NULL # to erase conflicting name
           ans$counts[2] <- ans$iter
           ans$fn.reduction <- NULL # so it does not interfere
           ans$iter<-NULL
           ans$gradient<-NULL # loss of information
        } else { # spg failed
		if (control$trace > 0) cat("spg failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        ## return(ans)
      }  # end if using spg
## --------------------------------------------
      else if (method == "ucminf") {
        ## Use ucminf routine
        if (is.null(control$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- control$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
        errmsg <- NULL
        if (control$have.bounds) {
              if (control$trace > 0) cat("ucminf cannot handle bounds\n")
              errmsg <- "ucminf cannot handle bounds\n"
            ##  stop("ucminf tried with bounds")
            ans <- list()
            class(ans)[1] <- "try-error"
        } else {
          
         uhessian <- 0 # Ensure hessian NOT computed
         ans <- try(ucminf::ucminf(par=spar, fn=efn, gr=egr, 
                   hessian = uhessian,  control=mcontrol, ...))
        }
        if (class(ans)[1] != "try-error") {
# From ucminf documentation:  convergence = 1 Stopped by small gradient (grtol).
#                                           2 Stopped by small step (xtol).
#                                           3 Stopped by function evaluation limit (maxeval).
#                                           4 Stopped by zero step from line search
#                                           -2 Computation did not start: length(par) = 0.
#                                           -4 Computation did not start: stepmax is too small.
#                                           -5 Computation did not start: grtol or xtol <= 0.
#                                           -6 Computation did not start: maxeval <= 0.
#                                           -7 Computation did not start: given Hessian not pos. definite.
#                             message: String with reason of termination.
		if (ans$convergence == 1 || ans$convergence == 2 || ans$convergence == 4) {
         		ans$convergence <- 0
		} 
                ans$par <- ans$par*pscale
        	ans$counts[1] <- ans$info[4]
        	ans$counts[2] <- ans$info[4] # calls fn and gr together
        	ans$info <- NULL # to erase conflicting name
        	ans$nitns <- NULL
                ans$hessian <- NULL
                ans$invhessian.lt <- NULL
		if (control$trace > 0) cat("ucminf message:",ans$message,"\n")
        } else { # ucminf failed
		if (control$trace > 0) cat("ucminf failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- errmsg
                ans$hessian <- NULL
        }
        uhessian <- NULL
        ## return(ans)
      }  ## end if using ucminf
## --------------------------------------------
      else if (method == "Rtnmin") { # Use Rtnmin routines 
	if (control$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}
	ans<-list() # ans not yet defined, so set as list
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
        if (is.null(egr)) { ## fixed msg below (referred to lbfgs) 170214
            if (control$trace > 0) cat("Rtnmin MUST have gradient provided\n")
            errmsg <- "Rtnmin MUST have gradient provided"
            class(ans)[1] <- "try-error"            
        }
	if (control$have.bounds) {
   	   ans <- try(Rtnmin::tnbc(x=spar, fgfun=nlmfn, lower=slower,
                upper=supper, trace=mcontrol$trace, ...))
	} else {
   	   ans <- try(Rtnmin::tn(x=spar, fgfun=nlmfn, trace=mcontrol$trace, ...))
	}
        if (class(ans)[1] == "try-error") {
        	if (control$trace>0) cat("Rtnmin failed for current problem \n")
                ans$convergence <- 9999 # failed in run
	        ans$message <- "Rtnmin failed fo current problem"        
                if (is.null(egr)) {
                   ans$convergence <- 9998
                   ans$message <- errmsg
                   ans$value <- 1234567E20
                } 
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA 
                ans$hessian <- NULL
        } else {
                ans$par <- ans$xstar*pscale
                ans$xstar <- NULL
                ans$value <- as.numeric(ans$f)
                ans$f <- NULL
                ans$g <- NULL
		ans$convergence <- ans$ierror
                ans$ierror <- NULL
	        ans$counts[1] <- ans$nfngr
	        ans$counts[2] <- ans$nfngr
                ans$nfngr <- NULL
                ans$hessian <- NULL
	        ans$message <- NA
        }
        ## return(ans)
      }  ## end if using Rtnmin
## --------------------------------------------
      else if (method == "bobyqa") {# Use bobyqa routine from minqa package
  	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        myrhobeg <- min(supper - slower)/3 # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        ans <- try(minqa::bobyqa(par=spar, fn=efn, lower=slower,
                upper=supper, control=mcontrol,...))
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
#                if (ans$feval > mcontrol$maxfun) {
#			ans$convergence <- 1 # too many evaluations
#                }
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
                ans$feval <- NULL
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using bobyqa
## --------------------------------------------
      else if (method == "uobyqa") {# Use uobyqa routine from minqa package
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        if (control$have.bounds) {
            warning("Cannot use uobyqa with bounds")
		if (control$trace > 0) cat("Cannot use uobyqa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
                ## return(ans)
        }
        ans <- try(minqa::uobyqa(par=spar, fn=efn, control=mcontrol,...))
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
#                if (ans$feval > mcontrol$maxfun) {
#			ans$convergence <- 1 # too many evaluations
#                }
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
                ans$feval <- NULL
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("uobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using uobyqa
## --------------------------------------------
      else if (method == "newuoa") {# Use newuoa routine from minqa package
        if (control$trace > 1) cat("Trying newuoa\n")
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        if (control$have.bounds) {
            warning("Cannot use newuoa with bounds")
		if (control$trace > 0) cat("Cannot use newuoa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
                ## return(ans)
        }
        ans <- try(minqa::newuoa(par=spar, fn=efn, control=mcontrol,...))
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
#                if (ans$feval > mcontrol$maxfun) {
#			ans$convergence <- 1 # too many evaluations
#                }
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
                ans$feval <- NULL
	        ans$counts[2] <- NA
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- defctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using newuoa
## --------------------------------------------
      else if (method == "nmkb") {# Use nmkb routine from dfoptim package
        if (any(par == lower) || any(par==upper)) {
           if (control$trace>0) cat("nmkb cannot start if on any bound \n")
           warning("nmkb() cannot be started if any parameter on a bound")
           ans <- list() # ans not yet defined, so set as list
           ans$value <- defctrl$badval
           ans$par <- rep(NA,npar)
           ans$convcode <- 9999 # failed in run - ?? consider special code for nmkb on bounds
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
           ans$hessian <- NULL
        } else { # ok to proceed with nmkb()
        if (! is.null(control$maxit)) { 
	   mcontrol$maxfeval <- control$maxit
	} else {
	   mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?
	}
         if (control$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         if (control$have.bounds) {
            ans <- try(dfoptim::nmkb(par=spar, fn=efn, lower = slower, 
              upper = supper, control=mcontrol, ...))
         } else {# 170919 explicit package in call
            ans <- try(dfoptim::nmk(par=spar, fn=efn, control=mcontrol, ...))
         }
         if (control$trace > 1) {
            cat("Outputting ans for nmkb:\n")
            print(ans)
         }

         if (class(ans)[1] != "try-error") {
           ans$value <- as.numeric(ans$value)
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval <- NULL
           ans$counts[2] <- NA
      	   ans$nitns <- NA # not used
           # What about 'restarts' and 'message'??
           warning(ans$message,"  Restarts for stagnation =",ans$restarts)
           ans$restarts <- NULL
           ans$hessian <- NULL
         } else {
           if (control$trace>0) cat("nmkb failed for current problem \n")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- defctrl$badval
           ans$par <- rep(NA,npar)
           ans$counts[1] <- NA
           ans$counts[2] <- NA
           ans$convergence <- 9999 # failed in run
           ans$message<-"Failed"
           ans$hessian <- NULL
         }
       } # end of check for parameter on bound
       ## return(ans)
     }  ## end if using nmkb
## --------------------------------------------
      else if (method == "hjkb") {# Use hjkb routine from dfoptim package
         if (control$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
         mcontrol$maxfeval <- control$maxfeval
         if (control$have.bounds) {
            ans <- try(dfoptim::hjkb(par=spar, fn=efn, lower = slower, 
                upper = supper, control=mcontrol, ...))
         } else {
            ans <- try(dfoptim::hjk(par=spar, fn=efn, control=mcontrol, ...))
         }
         if (class(ans)[1] != "try-error") {
           ans$value <- as.numeric(ans$value)
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval <- NULL
           ans$counts[2] <- NA
      	   ans$nitns <- NULL # not used
           ans$restarts <- NULL
           ans$hessian <- NULL
           ans$nitns <- NULL # loss of information
         } else {
            if (control$trace>0) cat("hjkb failed for current problem \n")
            ans <- list(value=defctrl$badval, par=rep(NA,npar), message="Failed",
                convergence=9999)
            ans$counts[1]<- NA
            ans$counts[2]<- NA 
            ans$hessian <- NULL
         }
         ## return(ans)
      }  ## end if using hjkb
## --------------------------------------------
      else if (method == "lbfgsb3") {# Use 2011 L-BFGS-B wrapper
        if (control$trace > 1) cat("lbfgsb3\n")
        mcontrol$trace <- control$trace
        if (control$trace < 1) {mcontrol$iprint <- -1} else {mcontrol$iprint <- control$trace} 
        # ?? use maxfevals rather than maxit for lbfgsb3 ?
        if (control$trace > 0) cat("lbfgsb3:control$have.bounds =",control$have.bounds,"\n")
        if (control$have.bounds) { ## Note call uses prm not par
            slower <- lower/pscale
            supper <- upper/pscale
            ans <- try(lbfgsb3::lbfgsb3(prm=spar, fn=efn, gr=egr, lower = slower, 
                upper = supper, control=mcontrol, ...)) # explicit pkg in call 170919
        } else {
            ans <- try(lbfgsb3::lbfgsb3(prm=spar, fn=efn, gr=egr, control=mcontrol, ...))
        }
        if (class(ans)[1] != "try-error") {
 ## Need to check these carefully??
            ans$convergence <- 0
            ans$par <- ans$prm*pscale
            ans$prm <- NULL
            ans$value<-as.numeric(ans$f)
            ans$f <- NULL
            ans$counts[1] <- ans$info$isave[34]
            ans$counts[2] <- ans$counts[1]
            ans$info <- NULL ## Note -- throwing away a lot of information
            ans$g <- NULL ## perhaps keep -- but how??
            ans$hessian <- NULL
            ans$message <- NA
            ans$niter <- NULL # loss of information
         } else {
            if (control$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- defctrl$badval
            ans$par<-rep(NA,npar)
            ans$convergence<-9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         ## return(ans)
      }  ## end if using lbfgsb3
## --------------------------------------------
      else if (method == "lbfgs") {# Use unconstrained method from lbfgs package
        if (control$trace > 1) cat("lbfgs\n")
        if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
        if (control$trace > 1) cat("lbfgs:control$have.bounds =",control$have.bounds,"\n")
        ans <- list() # to define the answer object
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
##      cat("in lbfgs section, control$have.bounds=",control$have.bounds,"\n")
        if (control$have.bounds) {
              cat("control$have.bounds seems TRUE\n")
              if (control$trace > 0) cat("lbfgs::lbfgs cannot handle bounds\n")
              errmsg <- "lbfgs::lbfgs cannot handle bounds\n"
            ##  stop("lbfgs::lbfgs tried with bounds")
            class(ans)[1] <- "try-error"            
        }
        if (is.null(egr)) {
            if (control$trace > 0) cat("lbfgs::lbfgs MUST have gradient provided\n")
            errmsg <- "lbfgs::lbfgs MUST have gradient provided\n"
            class(ans)[1] <- "try-error"            
        }
        if (class(ans)[1] == "undefined"){
            dotstuff <- list(...)
	    # cat("dotstuff:\n")
#	    print(dotstuff)
	    dotstuff$pscale <- pscale
	    dotstuff$fnscale <- fnscale
	    eopt <- list2env(dotstuff) # put it in an environment
	    # print(ls(eopt))
            ans <- try(lbfgs::lbfgs(efn, egr, vars=spar, 
                    environment=eopt, invisible=invisible))
        }
#        cat("interim answer:")
#        print(ans)
        if (class(ans)[1] != "try-error") {
        ## Need to check these carefully??
            ans$par <- ans$par*pscale
            ans$value <- ans$value*fnscale
            ans$counts[1] <- NA # lbfgs seems to have no output like this
            ans$counts[2] <- NA
         } else {
            if (control$trace>0) cat("lbfgs failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- defctrl$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
         ## return(ans)
      }  ## end if using lbfgs
  ## --------------------------------------------
  else if (method == "subplex") {# Use unconstrained method from subplex package
    if (control$trace > 1) cat("subplex\n")
    if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
    if (control$trace > 1) cat("subplex:control$have.bounds =",control$have.bounds,"\n")
    ans <- list() # to define the answer object
    errmsg <- NA
    if (control$trace > 0) warning("subplex has no trace mechanism")
    class(ans)[1] <- "undefined" # initial setting
    if (control$have.bounds) {
      cat("control$have.bounds seems TRUE\n")
      if (control$trace > 0) cat("subplex::subplex cannot handle bounds\n")
      errmsg <- "subplex::subplex cannot handle bounds\n"
      ##  stop("subplex::lbfgs tried with bounds")
      class(ans)[1] <- "try-error"            
    }
    if (class(ans)[1] == "undefined"){
          # Do we want parscale or use spar??
          cat("maxit =",defctrl$maxfeval,"\n")
          ans <- try(subplex::subplex(par=spar, fn=efn, control=list(maxit=defctrl$maxfeval)))
    }
    #        cat("interim answer:")
    #        print(ans)
    if (class(ans)[1] != "try-error") {
      ## Need to check these carefully??
      ans$par <- ans$par*pscale
      ans$value <- ans$value*fnscale
      ans$counts[1] <- ans$count
      ans$counts[2] <- NA
      ans$count <- NULL
      ccode <- ans$convergence
      if (ccode == -2) {ans$convergence <- 20}
          else {if (ccode == 0) {ans$convergence <- 0 } # not needed
               else {if (ccode == -1) {ans$convergence <- 1} else {ans$convergence <- 9999}}# unknown 
           }
    } else {
      if (control$trace>0) cat("subplex failed for current problem \n")
      ans<-list() # ans not yet defined, so set as list
      ans$value <- defctrl$badval
      ans$par <- rep(NA,npar)
      ans$convergence <- 9999 # failed in run
      if (is.null(egr)) ans$convergence <- 9998 # no gradient
      ans$counts[1] <- NA
      ans$counts[1] <- NA
      ans$hessian <- NULL
      if (! is.na(errmsg)) ans$message <- errmsg
    }
    ## return(ans)
  }  ## end if using subplex
## --------------------------------------------
## END OF optimrx extra methods
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD:", method, sep='')
             stop(errmsg, call.=FALSE)
      }
# Exit from routine
      if (savehess) { # compute hessian
         if (is.null(orig.gr)) {
            hess <- hessian(orig.fn, ans$par, ...) # from numDeriv
         } else { 
            hess <- jacobian(orig.gr, ans$par, ...) # use Jacobian of gradient
         }
      } else { hess <- NULL } # to ensure it is defined
      ans$hessian <- hess
      ans # last statement of routine
} ## end of optimrx
