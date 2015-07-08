optim <- function(par, ufn, ugr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {

# ?? time not used in output -- make it an attribute of ans??
# The structure has   par, value, counts, convergence, message, hessian

# Run a single method
## ?? check have.bounds?? or do again??
## 131027 ?? needed or in setup
  npar <- length(par)
  if (length(lower) == 1) lower<-rep(lower,npar)
  if (length(upper) == 1) upper<-rep(upper,npar)
## end 131027
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
        time <- system.time(ans <- try(optim(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper, 
                method=method, hessian=hessian, control=mcontrol, ...), silent=TRUE))[1]
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
        time <- system.time(ans <- try(nlminb(start=par, objective=ufn, gradient=ugr, lower=lower, 
		upper=upper, control=mcontrol,  ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		# Translate output to common format and names
        	ans$value<-ans$objective
	        ans$objective<-NULL
	        ans$counts[1] <- ans$evaluations[1]
        	ans$counts[2] <- ans$evaluations[2]
		ans$evaluations<-NULL # cleanup
	        ans$iterations<-NULL
	} else { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		if (control$trace>0) cat("nlminb failure\n")
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
        }
        return(ans)
      }  ## end if using nlminb
## --------------------------------------------
      else if (method == "nlm") { # Use stats package nlm routine
        if (!is.null(ugr)) { # ?? can ugr be null?
           tufn <- function(par, ...){
               f <- ufn(par, ...)
               g <- ugr(par, ...)
               attr(f,"gradient") <- g
               attr(f,"hessian") <- NULL # ?? maybe change later
               f
           }
        } else { stop("optimz::optim -- nlm -- we do not allow ugr = NULL") }
	if (! is.null(control$maxit) ) {iterlim <- control$maxit }
        else { iterlim <- 100 }
	print.level <- 0 
        if (! is.null(control$trace) && (control$trace > 0) ) {print.level <- 2 } 
        time <- system.time(ans <- try(nlm(f=tufn, p=par, iterlim=iterlim, 
                    print.level=print.level, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		if (ans$code == 1 || ans$code == 2 || ans$code == 3) ans$convergence <- 0
		if (ans$code == 4) ans$convvergence <- 1
                if (ans$code == 5) ans$convergence <- 5
        	# Translate output to common format
		ans$value <- ans$minimum
		ans$minimum <- NULL
                ans$par <- ans$estimate
		ans$estimate <- NULL
#        	ans$nitns <- ans$iterations
        	ans$iterations <- NULL
        	ans$counts[2] <- ans$nitns
                ans$counts[1] <- NA
	} else {
		if (control$trace > 0) cat("nlm failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
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
        time <- system.time(ans <- try(spg(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper,  
		control=mcontrol, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") { 
           ans$counts[1] <- ans$feval
           ans$feval<-NULL # to erase conflicting name
           ans$counts[2] <- ans$iter
           ans$iter<-NULL
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
        }
        return(ans)
      }  # end if using spg
## --------------------------------------------
      else if (method == "ucminf") {
        ## Use ucminf routine
        if (is.null(control$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- control$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
        if (hessian) uhessian <- 1 else uhessian <- 0
        time <- system.time(ans <- try(ucminf(par=par, fn=ufn, gr=ugr, 
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
        	ans$counts[1] <- ans$info[4]
        	ans$counts[2] <- ans$info[4] # calls fn and gr together
        	ans$info <- NULL # to erase conflicting name
        	ans$nitns <- NULL
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
        }
        uhessian <- NULL
        return(ans)
      }  ## end if using ucminf
## --------------------------------------------
      else if (method == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
	bdmsk <- bmchk(par, lower=lower, upper=upper)$bdmsk ## 131027 removed
        mcontrol <- control
	if (control$have.bounds) {
#           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminb(par=par, fn=ufn, gr=ugr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
#           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminu(par=par, fn=ufn, gr=ugr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
	        ans$fevals<-ans$counts[1]
	        ans$gevals<-ans$counts[2]
                ans$counts<-NULL
	#	ans$value<-ans$value 
        } else {
		if (control$trace>0) cat("Rcgmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
        }
        return(ans)
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (method == "Rtnmin") { # Use Rtnmin routines (ignoring masks)
	if (control$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}
         tufn <- function(par, ...){
            f <- ufn(par, ...)
            g <- ugr(par, ...)
            attr(f,"gradient") <- g
            attr(f,"hessian") <- NULL # ?? maybe change later
            f
        }
	if (control$have.bounds) {
   	   time <- system.time(ans <- try(tnbc(x=par, fgfun=tufn, gr=ugr, lower=lower,
                upper=upper, trace=mcontrol$trace, ...), silent=TRUE))[1]
	} else {
   	   time <- system.time(ans <- try(tn(x=par, fgfun=tufn, 
		  trace=mcontrol$trace, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
                ans$par <- ans$xstar
                ans$xstar <- NULL
                ans$value <- ans$f
                ans$f <- NULL
		ans$convergence <- ans$ierror
	        ans$counts[1] <- ans$nfngr
	        ans$counts[2] <- ans$nfngr
                ans$nfngr <- NULL
        } else {
		if (control$trace>0) cat("Rtnmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
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
   	   time <- system.time(ans <- try(Rvmminb(par=par, fn=ufn, gr=ugr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
   	   time <- system.time(ans <- try(Rvmminu(par=par, fn=ufn, gr=ugr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] == "try-error") {
            if (control$trace>0) cat("Rvmmin failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
        }
        return(ans)
      }  ## end if using Rvmmin
## --------------------------------------------
      else if (method == "bobyqa") {# Use bobyqa routine from minqa package
  	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        time <- system.time(ans <- try(bobyqa(par=par, fn=ufn, lower=lower, upper=upper, 
		control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convergence <- 1 # too many evaluations
                }
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
		ans$value<-ans$fval 
	      	ans$fval <- NULL # not used
        } else {
		if (control$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
        }
        return(ans)
      }  ## end if using bobyqa
## --------------------------------------------
      else if (method == "uobyqa") {# Use uobyqa routine from minqa package
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace

        time <- system.time(ans <- try(uobyqa(par=par, fn=ufn, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convergence <- 1 # too many evaluations
                }
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
		ans$value<-ans$fval 
	      	ans$fval <- NULL # not used
        } else {
		if (control$trace>0) cat("uobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
        }
        return(ans)
      }  ## end if using uobyqa
## --------------------------------------------
      else if (method == "newuoa") {# Use newuoa routine from minqa package
        cat("Trying newuoa\n")
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace

        time <- system.time(ans <- try(newuoa(par=par, fn=ufn, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convergence <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convergence <- 1 # too many evaluations
                }
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
		ans$value<-ans$fval 
	      	ans$fval <- NULL # not used
        } else {
		if (control$trace>0) cat("newuoa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
        }
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
        } else { # ok to proceed with nmkb()
        if (! is.null(control$maxit)) { 
		mcontrol$maxfeval <- control$maxit
	} else {
		mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
	}
         if (control$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         if (control$have.bounds) {
            time <- system.time(ans <- try(nmkb(par=par, fn=ufn, lower = lower, 
              upper = upper, control=mcontrol, ...), silent=TRUE))[1]
         } else {
            time <- system.time(ans <- try(nmk(par=par, fn=ufn, 
              control=mcontrol, ...), silent=TRUE))[1]
         }
         if (class(ans)[1] != "try-error") {
           ans$convcode <- ans$convergence
           ans$convergence <- NULL
           ans$value <- as.numeric(ans$value)
           ans$fevals <- ans$feval
           ans$feval <- NULL
           ans$gevals <- NA
           if (ans$fevals > mcontrol$maxfeval) { ## Note maxfeval defined above.
             ans$convcode <- 1 # too many evaluations
           }
      	   ans$nitns <- NA # not used
           # What about 'restarts' and 'message'??
           warning(ans$message,"  Restarts for stagnation =",ans$restarts)
           ans$message <- NULL
           ans$restarts <- NULL
         } else {
           if (control$trace>0) cat("nmkb failed for current problem \n")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- control$badval
           ans$par <- rep(NA,npar)
           ans$convcode <- 9999 # failed in run
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
         }
       } # end of check for parameter on bound
     }  ## end if using nmkb
## --------------------------------------------
      else 
      if (method == "hjkb") {# Use hjkb routine from dfoptim package
         if (control$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
         mcontrol$maxfeval <- control$maxfeval
         if (have.bounds) {
            time <- system.time(ans <- try(hjkb(par=par, fn=ufn, lower = lower, 
                upper = upper, control=mcontrol, ...), silent=TRUE))[1]
         } else {
            time <- system.time(ans <- try(hjk(par=par, fn=ufn, 
                control=mcontrol, ...), silent=TRUE))[1]
         }
         if (class(ans)[1] != "try-error") {
            ans$convcode <- ans$convergence
            if (ans$convcode == 1) ans$convcode=9999
            ans$convergence <- NULL
            ans$value <- as.numeric(ans$value)
            ans$fevals <- ans$feval
            ans$feval <- NULL
            if (ans$fevals > mcontrol$maxfeval) {
		ans$convcode <- 1 # too many evaluations
            }
            ans$gevals <- NA
            ans$nitns <- ans$niter
            ans$niter <- NULL
         } else {
            if (control$trace>0) cat("hjkb failed for current problem \n")
            ans <- list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convcode <- 9999 # failed in run
            ans$gevals <- NA 
            ans$nitns <- NA
         }
      }  ## end if using hjkb
## --------------------------------------------
      else 
      if (method == "lbfgsb3") {# Use 2011 L-BFGS-B wrapper
        if (control$trace > 1) cat("lbfgsb3\n")
        mcontrol <- control
#        mcontrol$maxfeval <- control$maxfeval
        # ?? use maxfevals rather than maxit for lbfgsb3 ??
        if (have.bounds) { ## Note call uses prm not par
            time <- system.time(ans <- try(lbfgsb3(par=par, fn=ufn, gr=ugr, lower = lower, 
                upper = upper, control=mcontrol, ...), silent=TRUE))[1]
        } else {
            time <- system.time(ans <- try(lbfgsb3(par=par, fn=ufn, gr=ugr, lower = -Inf, 
                upper = Inf, control=mcontrol, ...), silent=TRUE))[1]
        }
        if (class(ans)[1] != "try-error") {
 ## Need to check these carefully??
            ans$convcode <- 0
            ans$value<-as.numeric(ans$f)
            ans$fevals<-ans$info$isave[34]
            if (ans$fevals > mcontrol$maxfeval) {
               ans$convcode <- 1 # too many evaluations
            }
            ans$gevals<-ans$fevals
            ans$nitns<-ans$fevals
            ans$info <- NULL ##?? Note -- throwing away a lot of information
            ans$g <- NULL ##?? perhaps keep -- but how??
         } else {
            if (control$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par<-rep(NA,npar)
            ans$convcode<-9999 # failed in run
            ans$gevals<-NA 
            ans$nitns<-NA
         }
      }  ## end if using lbfgsb3
## --------------------------------------------
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD: ", method, sep='')
             stop(errmsg, call.=FALSE)
      }

} ## end of optimx.run
