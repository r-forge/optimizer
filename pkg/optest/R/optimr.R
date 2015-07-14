optimr <- function(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {

  if (is.null(control$trace)) control$trace <- 0
  npar <- length(par)
  if (is.null(control$parscale)) { pscale <- rep(1,npar) }
  else { pscale <- control$parscale }
  spar <- par/pscale # scaled parameters
  if (is.null(control$fnscale)) {
     fnscale <- 1
     if (! is.null(control$maximize) && control$maximize ) {fnscale <- -1}
  else if (! is.null(control$maximize)) {
          if ( (control$fnscale < 0) && maximize) {fnscale <- -1} # this is OK
          else stop("control$fnscale and control$maximize conflict")
       } # end ifelse
  } # end else

  efn <- function(spar, ...) {
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn(par, ...)
  }
  egr <- function(spar, ...) {
      par <- spar*pscale
      result <- gr(par, ...) * pscale
  }

  nlmfn <- function(spar, ...){
     f <- efn(spar, ...)
     g <- egr(spar, ...)
     attr(f,"gradient") <- g
     attr(f,"hessian") <- NULL # ?? maybe change later
     f
  }
# ?? do we want ehess ??    Not at 150714

# replacement for optim to minimize using a single method

# ?? time not used in output -- make it an attribute of ans??
# The structure has   par, value, counts, convergence, message, hessian

#  cat("optimr: control$trace =",control$trace,"\n")
#  print(par)
#  cat("fval from efn:")
#  print(efn(spar,pscale, fnscale, ...))
#  if (! is.null(gr)) {
#  cat("gval from efn:")
#  print(egr(par, pscale, fnscale, ...))
#  }
#  print(method)
# tmp <- readline("on to the run stage of optimr")

# Run a single method
  if (is.null(control$have.bounds)) {
     have.bounds <- bmchk(par, lower, upper)$bounds
     control$have.bounds <- have.bounds # for safety
     cat("optimr -- have.bounds =",have.bounds,"\n")
  } else have.bounds <- control$have.bounds

# expand bounds
  if (length(lower) == 1) lower<-rep(lower,npar)
  if (length(upper) == 1) upper<-rep(upper,npar)

  mcontrol <- list() # define the control list

# Methods from optim()
      if (method=="Nelder-Mead" || 
          method == "BFGS" || 
          method == "L-BFGS-B" || 
          method == "CG" || 
          method == "SANN") {
        # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
        mcontrol$maxit <- control$maxit 
        mcontrol$trace <- control$trace
	mcontrol$parscale <- control$parscale # Use internal scaling
# Note: hessian always FALSE

        if (have.bounds) {
        time <- system.time(ans <- try(optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, 
                method=method, hessian=FALSE, control=mcontrol, ...), silent=TRUE))[1]
        } else {
        time <- system.time(ans <- try(optim(par=par, fn=fn, gr=gr, 
                method=method, hessian=FALSE, control=mcontrol, ...), silent=TRUE))[1]
        }

        # The time is the index=1 element of the system.time for the process, 
        # which is a 'try()' of the regular optim() function
        if (class(ans)[1] == "try-error") { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		if (control$trace>0) cat("optim method failure\n")
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
        } # otherwise ans is OK and we return it
        return(ans) # to ensure we return
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
        time <- system.time(ans <- try(nlminb(start=spar, objective=efn, gradient=egr, lower=lower, 
		upper=upper, control=mcontrol,  ...), silent=TRUE))[1]
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        return(ans)
      }  ## end if using nlminb
## --------------------------------------------
      else if (method == "nlm") { # Use stats package nlm routine
        if (is.null(gr)) { stop("optimr -- nlm -- we do not allow gr = NULL") }
        cat("test nlmfn:")
        ffval <- nlmfn(spar)
        print(ffval)
        tmp <- readline("cont.")
	if (! is.null(control$maxit) ) {iterlim <- control$maxit }
        else { iterlim <- 100 }
	print.level <- 0 
        if (! is.null(control$trace) && (control$trace > 0) ) {print.level <- 2 } 
        time <- system.time(ans <- try(nlm(f=nlmfn, p=spar, iterlim=iterlim, 
                    print.level=print.level, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		if (ans$code == 1 || ans$code == 2 || ans$code == 3) ans$convergence <- 0
		if (ans$code == 4) ans$convvergence <- 1
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
	}
        print.level <- NULL # clean up
        return(ans)
      } # end if using nlm
## --------------------------------------------
      else if (method == "spg") { # Use BB package routine spg as minimizer
        mcontrol$maximize <- NULL # Use external maximization approach
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
        if (control$trace > 0) { 
            mcontrol$trace <- TRUE
            if (control$trace > 1) mcontrol$triter <- 1 # default is 10
        } else { mcontrol$trace <- FALSE }
        time <- system.time(ans <- try(spg(par=spar, fn=efn, gr=egr, lower=lower, upper=upper,  
		control=mcontrol, ...), silent=TRUE))[1]
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
		if (control$trace > 0) cat("nlm failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        return(ans)
      }  # end if using spg
## --------------------------------------------
      else if (method == "ucminf") {
        ## Use ucminf routine
        if (is.null(control$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- control$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
#        if (hessian) uhessian <- 1 else uhessian <- 0
         uhessian <- 0 # Ensure hessian NOT computed
        time <- system.time(ans <- try(ucminf(par=spar, fn=efn, gr=egr, 
                 hessian = uhessian,  control=mcontrol, ...), silent=TRUE))[1]
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        uhessian <- NULL
        return(ans)
      }  ## end if using ucminf
## --------------------------------------------
      else if (method == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
	bdmsk <- bmchk(par, lower=lower, upper=upper)$bdmsk ## 131027 removed
        mcontrol$trace <- control$trace
#?? others
	if (control$have.bounds) {
#           if (is.null(gr)) gr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminb(par=spar, fn=efn, gr=egr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
#           if (is.null(gr)) gr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminu(par=spar, fn=efn, gr=egr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (control$trace>0) cat("Rcgmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        return(ans)
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (method == "Rtnmin") { # Use Rtnmin routines (ignoring masks)
	if (control$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}

	if (control$have.bounds) {
   	   time <- system.time(ans <- try(tnbc(x=spar, fgfun=nlmfn, gr=gr, lower=lower,
                upper=upper, trace=mcontrol$trace, ...), silent=TRUE))[1]
	} else {
   	   time <- system.time(ans <- try(tn(x=spar, fgfun=nlmfn, 
		  trace=mcontrol$trace, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
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
        } else {
		if (control$trace>0) cat("Rtnmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        tufn <- NULL # 140902 -- clear function
        return(ans)
      }  ## end if using Rtnmin
## --------------------------------------------
      else if (method == "Rvmmin") { # Use Rvmmin routine (ignoring masks??)
	bdmsk<-bmchk(par, lower=lower, upper=upper)$bdmsk ## 131027 removed
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
	mcontrol$trace <- control$trace # 140902 Note no check on validity of values
	if (control$have.bounds) {
   	   time <- system.time(ans <- try(Rvmminb(par=par, fn=fn, gr=gr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
   	   time <- system.time(ans <- try(Rvmminu(par=par, fn=fn, gr=gr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
            ans$par <- ans$par*pscale
            ans$bdmsk <- NULL
        } else {
            if (control$trace>0) cat("Rvmmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        return(ans)
      }  ## end if using Rvmmin
## --------------------------------------------
      else if (method == "bobyqa") {# Use bobyqa routine from minqa package
  	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        time <- system.time(ans <- try(bobyqa(par=par, fn=fn, lower=lower, upper=upper, 
		control=mcontrol,...), silent=TRUE))[1]
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        return(ans)
      }  ## end if using bobyqa
## --------------------------------------------
      else if (method == "uobyqa") {# Use uobyqa routine from minqa package
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace

        time <- system.time(ans <- try(uobyqa(par=par, fn=fn, control=mcontrol,...), silent=TRUE))[1]
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        return(ans)
      }  ## end if using uobyqa
## --------------------------------------------
      else if (method == "newuoa") {# Use newuoa routine from minqa package
        cat("Trying newuoa\n")
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace

        time <- system.time(ans <- try(newuoa(par=par, fn=fn, control=mcontrol,...), silent=TRUE))[1]
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
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        return(ans)
      }  ## end if using newuoa
## --------------------------------------------
      else 
      if (method == "nmkb") {# Use nmkb routine from dfoptim package
        if (any(par == lower) || any(par==upper)) {
           if (control$trace>0) cat("nmkb cannot start if on any bound \n")
           warning("nmkb() cannot be started if any parameter on a bound")
           ans <- list() # ans not yet defined, so set as list
           ans$value <- control$badval
           ans$par <- rep(NA,npar)
           ans$convcode <- 9999 # failed in run - ?? need special code for nmkb on bounds
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
           ans$hessian <- NULL
        } else { # ok to proceed with nmkb()
        if (! is.null(control$maxit)) { 
		mcontrol$maxfeval <- control$maxit
	} else {
		mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
	}
         if (control$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         if (control$have.bounds) {
            time <- system.time(ans <- try(nmkb(par=spar, fn=efn, lower = lower, 
              upper = upper, control=mcontrol, ...), silent=TRUE))[1]
         } else {
            time <- system.time(ans <- try(nmk(par=spar, fn=efn, 
              control=mcontrol, ...), silent=TRUE))[1]
         }
         cat("Outputting ans for nmkb:\n")
         print(ans)

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
           ans$value <- control$badval
           ans$par <- rep(NA,npar)
           ans$convergence <- 9999 # failed in run
           ans$message<-"Failed"
           ans$hessian <- NULL
         }
       } # end of check for parameter on bound
       return(ans)
     }  ## end if using nmkb
## --------------------------------------------
      else 
      if (method == "hjkb") {# Use hjkb routine from dfoptim package
         if (control$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
         mcontrol$maxfeval <- control$maxfeval
         if (have.bounds) {
            time <- system.time(ans <- try(hjkb(par=spar, fn=efn, lower = lower, 
                upper = upper, control=mcontrol, ...), silent=TRUE))[1]
         } else {
            time <- system.time(ans <- try(hjk(par=spar, fn=efn, 
                control=mcontrol, ...), silent=TRUE))[1]
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
            ans <- list(value=control$badval, par=rep(NA,npar), message="Failed",
                convergence=9999)
            ans$counts[1]<- NA
            ans$counts[2]<- NA 
            ans$hessian <- NULL
         }
         return(ans)
      }  ## end if using hjkb
## --------------------------------------------
      else 
      if (method == "lbfgsb3") {# Use 2011 L-BFGS-B wrapper
        if (control$trace > 1) cat("lbfgsb3\n")
        mcontrol$trace <- control$trace
        mcontrol$iprint <- 1 # ???
#        mcontrol$maxfeval <- control$maxfeval
        # ?? use maxfevals rather than maxit for lbfgsb3 ??
        cat("have.bounds =",have.bounds,"\n")
        if (have.bounds) { ## Note call uses prm not par
            time <- system.time(ans <- try(lbfgsb3(prm=spar, fn=efn, gr=egr, lower = lower, 
                upper = upper, control=mcontrol, ...), silent=TRUE))[1]
        } else {
            time <- system.time(ans <- try(lbfgsb3(prm=spar, fn=efn, gr=egr,  
                control=mcontrol, ...), silent=TRUE))[1]
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
            ans$info <- NULL ##?? Note -- throwing away a lot of information
            ans$g <- NULL ##?? perhaps keep -- but how??
            ans$hessian <- NULL
            ans$message <- NA
            ans$niter <- NULL # loss of information
         } else {
            if (control$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par<-rep(NA,npar)
            ans$convergence<-9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         return(ans)
      }  ## end if using lbfgsb3
## --------------------------------------------
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD: ", method, sep='')
             stop(errmsg, call.=FALSE)
      }

} ## end of optimx.run
