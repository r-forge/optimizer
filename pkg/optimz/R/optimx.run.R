optimx.run <- function(par, ufn, ugr=NULL, uhess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, itnmax=NULL, hessian=FALSE,
            ctrl, ...) {
# ?? change all defaults to be brought in from outside -- should have method=NULL here

cat("150708 - check ctrl\n")
print(ctrl)

# Run methods
  have.bounds<-ctrl$have.bounds
  ctrl$have.bounds<-NULL ## or we get errors in optim()
  npar<-length(par)
## 131027 -- ?? do we need this or can we move to setup or check?
  if (length(lower) == 1) lower<-rep(lower,npar)
  if (length(upper) == 1) upper<-rep(upper,npar)
## end 131027
  nmeth<-length(method)
  pstring<-names(par)
  if (is.null(pstring)) {
    pstring <- NULL
    for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
  }  
  cnames <- c(pstring, "value", "fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes")
  ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+8)
  colnames(ans.ret)<-cnames
  row.names(ans.ret)<-method
  ans.details <- list()
  ansout <- NULL # ensure NULL if we have no parameters or no successes
  numgrad<-FALSE
  # 130924 -- need to fix -- some methods can handle null gradient
  for (i in 1:nmeth) { # loop over the methods
      meth <- method[i] # extract the method name
      conv <- -1 # indicate that we have not yet converged
      # 20100608 - take care of polyalgorithms
      if (! is.null(itnmax) ) {
	if (length(itnmax) == 1) {ctrl$maxit <- itnmax} # Note we will execute this FIRST
        else {if (length(itnmax) != nmeth) { 
		stop("Length of itnmax =",length(itnmax)," but should be ",nmeth) }
              else { ctrl$maxit<-itnmax[i] }
        }
        if (ctrl$follow.on && (ctrl$trace>0)) cat("Do ",ctrl$maxit," steps of ")
      }
      if (ctrl$trace>0) cat("Method: ", meth, "\n") # display the method being used
      # Extract control information e.g., trace
      # 20100215: Note that maxit needs to be defined other than 0 e.g., for ucminf
      # create local control list for each single method -- this is one of the key issues for optimx
      # 20140902 big change for mcontrol -- start NULL and then sort out for each method
      mcontrol <- NULL
      # not used in any methods -- it is here for the scale check of parameters and bounds above
      ans.ret[i, "value"] <- .Machine$double.xmax # to ensure defined for sort
# Methods from optim()
      if (meth == "SANN") stop("SANN is not permitted within optimx") # 20140902
      if (meth=="Nelder-Mead" || meth == "BFGS" || meth == "L-BFGS-B" || meth == "CG") {
        # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
        mcontrol$maxit <- ctrl$maxit 
        mcontrol$trace <- ctrl$trace # for optim() at least
        ## add explicit methopts ??
        time <- system.time(ans <- try(optim(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper, 
                method=meth, control=mcontrol, ...), silent=TRUE))[1]
        # The time is the index=1 element of the system.time for the process, 
        # which is a 'try()' of the regular optim() function
        if (class(ans)[1] != "try-error") { 
                ans$convcode <- ans$convergence
                ans$convergence <- NULL
		#      convergence: An integer code. '0' indicates successful convergence.
#                if (meth=="SANN") ans$convcode = 1 # always the case for SANN (but it reports 0!)
	        ans$fevals<-ans$counts[1] # save function and gradient count information
	        ans$gevals<-ans$counts[2]
        	ans$counts<-NULL # and erase the counts element now data is saved
	} else { # bad result -- What to do?
		ans<-list(fevals=NA) # ans not yet defined, so set as list
                ans$convcode<-9999 # failed in run
		if (ctrl$trace>0) cat("optim function evaluation failure\n")
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
	        ans$fevals<-NA # save function and gradient count information
	        ans$gevals<-NA
        }
      	ans$nitns<-NA # not used
      }   # end if using optim() methods
## --------------------------------------------
      else if (meth == "nlminb") {
        # Here we use portLib routine nlminb rather than optim as our minimizer
        mcontrol$iter.max<-mcontrol$maxit # different name for iteration limit in this routine
        mcontrol$maxit<-NULL
        mcontrol$abs.tol<-0 # To fix issues when minimum is less than 0. 20100711
	if ( is.null(mcontrol$trace) || is.na(mcontrol$trace) || mcontrol$trace == 0) { 
		mcontrol$trace = 0
	} else { 
		mcontrol$trace = 1 # this is EVERY iteration. nlminb trace is freq of reporting.
	}
        time <- system.time(ans <- try(nlminb(start=par, objective=ufn, gradient=ugr, lower=lower, 
		upper=upper, control=mcontrol,  ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- ans$convergence
	        # Translate output to common format and names
        	ans$value<-ans$objective
	        ans$objective<-NULL
	        ans$fevals<-ans$evaluations[1]
        	ans$gevals<-ans$evaluations[2]
		ans$evaluations<-NULL # cleanup
        	ans$nitns<-ans$iterations
	        ans$iterations<-NULL
	} else { # bad result -- What to do?
		ans<-list(fevals=NA) # ans not yet defined, so set as list
                ans$convcode<-9999 # failed in run
		if (ctrl$trace>0) cat("nlminb function evaluation failure\n")
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
        	ans$nitns<-NA # not used
                ans$gevals<-NA ## ?? missing 130929
        }
        ans$convergence<-NULL
      }  ## end if using nlminb
## --------------------------------------------
      else if (meth == "nlm") { # Use stats package nlm routine
#        tufn <- ufn # we don't want to change the user function object, so create a copy 
##?? This is giving attribute to the function rather than its value
#        if (!is.null(ugr)) {
#	   attr(tufn, "gradient") <- ugr(par, ...) 
#           if (!is.null(uhess)) {
##		if (ctrl$maximize) attr(tufn, "hessian") <- - uhess(par, ...) 
#                     # Only change attibute if gr defined too
##			else attr(tufn, "hessian") <-  uhess(par, ...)
#		attr(tufn, "hessian") <-  uhess(par, ...)
#           } else attr(tufn, "hessian") <- NULL
#        }
# 140902 try again
        if (!is.null(ugr)) { # ?? can ugr be null?
           tufn <- function(par, ...){
               f <- ufn(par, ...)
               g <- ugr(par, ...)
               attr(f,"gradient") <- g
               attr(f,"hessian") <- NULL # ?? maybe change later
               f
           }
        } else {
           tufn <- function(par, ...){
               f <- ufn(par, ...)
               g <- grfwd(ufn, par, ...) ##?? should allow choice
               attr(f,"gradient") <- g
               attr(f,"hessian") <- NULL # ?? maybe change later
               f
           }
        }
	## 091215 added control for iteration limit
	if (! is.null(ctrl$maxit)) { iterlim <- ctrl$maxit } 
	else { iterlim <- 100 }
	print.level <- 0 
        if (! is.null(ctrl$trace) && (ctrl$trace > 0) ) {print.level <- 2 } 
        # 20140902 Note: May want to check validity of values??
	# 110121 -- need to put tufn NOT ufn in call 
        # 140902 Note that nlm has arguments (not control arguments) such as gradtol. That needs special
	# attention if we decide to include??
    
        time <- system.time(ans <- try(nlm(f=tufn, p=par, iterlim=iterlim, 
                    print.level=print.level, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- ans$code
		if (ans$convcode == 1 || ans$convcode == 2 || ans$convcode == 3) ans$convcode <- 0
		if (ans$convcode == 4) ans$convcode <- 1
        	# Translate output to common format
		ans$value<-ans$minimum
                ans$par<-ans$estimate
		ans$estimate<-NULL
		ans$minimum<-NULL
        	ans$nitns<-ans$iterations
        	ans$fevals<-ans$nitns
        	ans$gevals<-ans$nitns # Artificial copy of these values
        	ans$iterations<-NULL
	} else {
		if (ctrl$trace > 0) cat("nlm failed for this problem\n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode<-9999 # failed in run
        	ans$gevals<-NA 
        	ans$nitns<-NA
	}
        print.level <- NULL # clean up
      } # end if using nlm
## --------------------------------------------
      else if (meth == "spg") { # Use BB package routine spg as minimizer
        mcontrol$maximize <- NULL # Use external maximization approach
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval
        if (ctrl$trace > 0) { 
            mcontrol$trace <- TRUE
            if (ctrl$trace > 1) mcontrol$triter <- 1 # default is 10
        } else { mcontrol$trace <- FALSE }
        time <- system.time(ans <- try(spg(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper,  
		control=mcontrol, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") { 
   	   ans$convcode <- ans$convergence
           ans$fevals<-ans$feval
           ans$feval<-NULL # to erase conflicting name
           ans$gevals<-NA # ??fixup needed
           ans$nitns<-ans$iter
           ans$iter<-NULL
        } else { # spg failed
		if (ctrl$trace > 0) cat("spg failed for this problem\n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode<-9999 # failed in run
        	ans$gevals<-NA 
        	ans$nitns<-NA
        }
        ans$convergence<-NULL
##	if (ctrl$maximize) {
##		ans$value <- -ans$value
##	}
      }  # end if using spg
## --------------------------------------------
      else if (meth == "ucminf") {
        ## Use ucminf routine
        if (is.null(ctrl$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- ctrl$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
        time <- system.time(ans <- try(ucminf(par=par, fn=ufn, gr=ugr,  control=mcontrol, ...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- ans$convergence
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
		if (ans$convcode == 1 || ans$convcode == 2 || ans$convcode == 4) {
         		ans$convcode <- 0
		} else {
			ans$convcode <- ans$convergence
		} # Termination criteria are tricky here!  How to determine successful convergence?
        	ans$fevals<-ans$info[4]
        	ans$gevals<-ans$info[4] # calls fn and gr together
        	ans$info<-NULL # to erase conflicting name
        	ans$nitns<-NA
		if (ctrl$trace > 0) cat("ucminf message:",ans$message,"\n")
        } else { # ucminf failed
		if (ctrl$trace > 0) cat("ucminf failed for this problem\n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode<-9999 # failed in run
        	ans$gevals<-NA 
        	ans$nitns<-NA
        }
        ans$convergence<-NULL
      }  ## end if using ucminf
## --------------------------------------------
      else if (meth == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
	bdmsk <- bmchk(par, lower=lower, upper=upper)$bdmsk ## 131027 removed
        mcontrol <- ctrl
#	mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values
#        mcontrol$maxit <- ctrl$maxit
	if (have.bounds) {
           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminb(par=par, fn=ufn, gr=ugr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rcgminu(par=par, fn=ufn, gr=ugr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
		ans$convcode <- ans$convergence
	        ans$fevals<-ans$counts[1]
	        ans$gevals<-ans$counts[2]
                ans$counts<-NULL
	#	ans$value<-ans$value 
        } else {
		if (ctrl$trace>0) cat("Rcgmin failed for current problem \n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode<-9999 # failed in run
                ans$fevals<-NA
        	ans$gevals<-NA 
        }
       	ans$nitns<-NA
        ans$convergence<-NULL
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (meth == "Rtnmin") { # Use Rtnmin routines (ignoring masks)
	if (ctrl$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}
        if (!is.null(ugr)) { # ?? can ugr be null?
           tufn <- function(par, ...){
               f <- ufn(par, ...)
               g <- ugr(par, ...)
               attr(f,"gradient") <- g
##141010               attr(f,"hessian") <- NULL # ?? maybe change later
              f
           }
        } else {
           stop("no gradient for tn")
##141010           tufn <- function(par, ...){
##141010               f <- ufn(par, ...)
##141010               g <- grfwd(ufn, par, ...) ##?? should allow choice
##141010               attr(f,"gradient") <- g
##141010               attr(f,"hessian") <- NULL # ?? maybe change later
##141010               f
##141010           }
        }
	if (have.bounds) {
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
		ans$convcode <- ans$ierror
	        ans$fevals <- ans$nfngr
	        ans$gevals <- ans$nfngr
                ans$nfngr <- NULL
        } else {
		if (ctrl$trace>0) cat("Rtnmin failed for current problem \n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par <- rep(NA,npar)
                ans$convcode <- 9999 # failed in run
                ans$fevals <- NA
        	ans$gevals <- NA 
        }
       	ans$nitns <- NA
        ans$convergence <- NULL
        tufn <- NULL # 140902 -- clear function
      }  ## end if using Rtnmin
## --------------------------------------------
      else if (meth == "Rvmmin") { # Use Rvmmin routine (ignoring masks??)
	bdmsk<-bmchk(par, lower=lower, upper=upper)$bdmsk ## 131027 removed
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval
	mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values

	if (have.bounds) {
           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rvmminb(par=par, fn=ufn, gr=ugr, lower=lower,
                upper=upper, bdmsk=bdmsk, control=mcontrol, ...), silent=TRUE))[1]
	} else {
           if (is.null(ugr)) ugr<-"grfwd" ##JN
   	   time <- system.time(ans <- try(Rvmminu(par=par, fn=ufn, gr=ugr, 
		control=mcontrol, ...), silent=TRUE))[1]
	}
        if (class(ans)[1] != "try-error") {
		ans$convcode <- ans$convergence
	        ans$fevals<-ans$counts[1]
	        ans$gevals<-ans$counts[2]
        } else {
		if (ctrl$trace>0) cat("Rvmmin failed for current problem \n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode <- 9999 # failed in run
        	ans$gevals <- NA 
        }
       	ans$nitns <- NA # not used
        ans$convergence <- NULL
      }  ## end if using Rvmmin
## --------------------------------------------
      else if (meth == "bobyqa") {# Use bobyqa routine from minqa package
  	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace
        time <- system.time(ans <- try(bobyqa(par=par, fn=ufn, lower=lower, upper=upper, 
		control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convcode <- 1 # too many evaluations
                }
	        ans$fevals<-ans$feval
	        ans$gevals<-NA
		ans$value<-ans$fval 
	      	ans$nitns<-NA # not used
        } else {
		if (ctrl$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
                ans$convcode<-9999 # failed in run
        	ans$gevals<-NA 
        	ans$nitns<-NA
        }
      }  ## end if using bobyqa
## --------------------------------------------
      else if (meth == "uobyqa") {# Use uobyqa routine from minqa package
	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace

        time <- system.time(ans <- try(uobyqa(par=par, fn=ufn, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convcode <- 1 # too many evaluations
                }
	        ans$fevals <- ans$feval
	        ans$gevals <- NA
		ans$value <- ans$fval 
	      	ans$nitns <- NA # not used
        } else {
		if (ctrl$trace>0) cat("uobyqa failed for current problem \n")
		ans <- list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par <- rep(NA,npar)
                ans$convcode <- 9999 # failed in run
        	ans$gevals <- NA 
        	ans$nitns <- NA
        }
      }  ## end if using uobyqa
## --------------------------------------------
      else if (meth == "newuoa") {# Use newuoa routine from minqa package
        cat("Trying newuoa\n")
	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace

        time <- system.time(ans <- try(newuoa(par=par, fn=ufn, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans)[1] != "try-error") {
		ans$convcode <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$convcode <- 1 # too many evaluations
                }
	        ans$fevals <- ans$feval
	        ans$gevals <- NA
		ans$value <- ans$fval 
	      	ans$nitns <- NA # not used
        } else {
		if (ctrl$trace>0) cat("newuoa failed for current problem \n")
		ans <- list(fevals=NA) # ans not yet defined, so set as list
		ans$value <- ctrl$badval
		ans$par <- rep(NA,npar)
                ans$convcode <- 9999 # failed in run
        	ans$gevals <- NA 
        	ans$nitns <- NA
        }
      }  ## end if using newuoa
## --------------------------------------------
      else 
      if (meth == "nmkb") {# Use nmkb routine from dfoptim package
        if (any(par == lower) || any(par==upper)) {
           if (ctrl$trace>0) cat("nmkb cannot start if on any bound \n")
           warning("nmkb() cannot be started if any parameter on a bound")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- ctrl$badval
           ans$par <- rep(NA,npar)
           ans$convcode <- 9999 # failed in run - ?? need special code for nmkb on bounds
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
        } else { # ok to proceed with nmkb()
        if (! is.null(ctrl$maxit)) { 
		mcontrol$maxfeval <- ctrl$maxit
	} else {
		mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
	}
         if (ctrl$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         if (have.bounds) {
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
           if (ctrl$trace>0) cat("nmkb failed for current problem \n")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- ctrl$badval
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
      if (meth == "hjkb") {# Use hjkb routine from dfoptim package
         if (ctrl$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
#?? can't be null         if (! is.null(ctrl$maxfit)) { 
         mcontrol$maxfeval <- ctrl$maxfeval
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
            if (ctrl$trace>0) cat("hjkb failed for current problem \n")
            ans <- list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- ctrl$badval
            ans$par <- rep(NA,npar)
            ans$convcode <- 9999 # failed in run
            ans$gevals <- NA 
            ans$nitns <- NA
         }
      }  ## end if using hjkb
## --------------------------------------------
      else 
      if (meth == "lbfgsb3") {# Use 2011 L-BFGS-B wrapper
        if (ctrl$trace > 1) cat("lbfgsb3\n")
        mcontrol <- ctrl
#        mcontrol$maxfeval <- ctrl$maxfeval
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
            if (ctrl$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- ctrl$badval
            ans$par<-rep(NA,npar)
            ans$convcode<-9999 # failed in run
            ans$gevals<-NA 
            ans$nitns<-NA
         }
      }  ## end if using lbfgsb3
## --------------------------------------------
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD: ", meth, sep='')
             stop(errmsg, call.=FALSE)
      }
## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      if (ctrl$trace>0) { cat("Post processing for method ",meth,"\n") }
      if (exists("ans$message")) {
           amsg<-ans$message
           ans$message <- NULL
      } else { amsg <- "none" }
      ngatend <- NA
      nhatend <- NA
      hev <- NA
      if ( ctrl$save.failures || (ans$convcode == 0) ){# Save soln if converged or directed to save
          if (ctrl$trace && ans$convcode==0) cat("Successful convergence! \n") 
          # Testing final soln. Use numDeriv for gradient & Hessian; compute Hessian eigenvalues
          hessOK<-FALSE # indicator for later
          gradOK<-FALSE
          if ((ctrl$kkt || hessian) && (ans$convcode != 9999)) {
              if (ctrl$trace>0) cat("Compute Hessian approximation at finish of ",method[i],"\n")
              if (!is.null(uhess)){ # check if we have analytic hessian 
                 nhatend<-try(uhess(ans$par, ...), silent=TRUE)
                 if (class(nhatend) != "try-error") {
                    hessOK<-TRUE
                 }
              } else {
                 if (is.null(ugr)) {
                     nhatend<-try(hessian(ufn, ans$par, ...), silent=TRUE) # change 20100711
                 } else {
                     nhatend<-try(jacobian(ugr,ans$par, ...), silent=TRUE) # change 20100711
                 } # numerical hessian at "solution"
                 if (class(nhatend) != "try-error") {
                    hessOK<-TRUE
                 }
              } # end hessian calculation
          } # end test if hessian computed
          ans$kkt1<-NA
          ans$kkt2<-NA
          if ((hessian || ctrl$kkt) && (ans$convcode != 9999)) {# avoid test when method failed
             if (ctrl$trace>0) cat("Compute gradient approximation at finish of ",method[i],"\n")
             if (is.null(ugr)) {
                 ngatend<-try(grad(ufn, ans$par, ...), silent=TRUE) # change 20100711
             } else {
                 ngatend<-try(ugr(ans$par, ...), silent=TRUE) # Gradient at solution # change 20100711
             }
             if (class(ngatend) != "try-error") gradOK<-TRUE # 100215 had == rather than != here
             if ( (! gradOK) && (ctrl$trace>0)) cat("Gradient computation failure!\n") 
             if (gradOK) {
                # test gradient
                ans$kkt1<-(max(abs(ngatend)) <= ctrl$kkttol*(1.0+abs(ans$value)) ) # ?? sensible?
                if (hessOK) {
                   # For bounds constraints, we need to "project" the gradient and Hessian
                   bset<-sort(unique(c(which(ans$par<=lower), which(ans$par>=upper))))
                   nbds<-length(bset) # number of elements nulled by bounds
                   # Note: we assume that we are ON, not over boundary, 
                   # but use <= and >=. No tolerance is used.
                   ngatend[bset]<-0 # "project" the gradient
                   nhatend[bset,] <-0
                   nhatend[, bset] <-0 # and the Hessian
                   if (!isSymmetric(nhatend, tol=ctrl$hessasymtol)) {
                      # hessOK<-FALSE ## Assume we will keep hessian after symmetrizing
                      asym <- sum(abs(t(nhatend) - nhatend))/sum(abs(nhatend))
                      asw <- paste("Hessian is reported non-symmetric with asymmetry ratio ", 
                      asym, sep = "")
                      if (ctrl$trace > 1) cat(asw, "\n")
                      if (ctrl$dowarn) warning(asw)
                      if (ctrl$trace > 1) cat("Force Hessian symmetric\n")
                      if (ctrl$dowarn) warning("Hessian forced symmetric", call. = FALSE)
                      nhatend <- 0.5 * (t(nhatend) + nhatend)
                   }  # end symmetry test
                   hev<- try(eigen(nhatend)$values, silent=TRUE) # 091215 use try in case of trouble
                   if (ctrl$kkt){
   	              if (class(hev) != "try-error") {# answers are OK, check Hessian properties
                         if (any(is.complex(hev))){
                            hessOK<-FALSE
                            cat("Complex eigenvalues found for method =",meth,"\n")
                            cat("coefficients for function value", ans$value," :\n")
                            print(ans$par)
                            dput(nhatend, file="badhess.txt")
                            warning("Complex eigenvalues found for method =",meth)
                         }
                         if (hessOK) {
                            negeig<-(hev[npar] <= (-1)*ctrl$kkt2tol*(1.0+abs(ans$value))) 
                            evratio<-hev[npar-nbds]/hev[1]
                            # If non-positive definite, then there have zero evs (from the projection)
                            # in the place of the "last" eigenvalue and we'll have singularity.
                            # WARNING: Could have a weak minimum if semi-definite.
                            ans$kkt2<- (evratio > ctrl$kkt2tol) && (! negeig)
                         }
                      } else {
                         warnstr<-paste("Eigenvalue failure after method ",method[i],sep='')
                         if (ctrl$dowarn) warning(warnstr)
                         if (ctrl$trace>0) {
                            cat("Hessian eigenvalue calculation failure!\n")
                            print(nhatend)
                         }
                      }
                   } # kkt2 evaluation
                } else { # computing Hessian has failed
                   warnstr<-paste("Hessian not computable after method ",method[i],sep='')
                   if (ctrl$dowarn) warning(warnstr)
                   if (ctrl$trace>0) cat(warnstr,"\n") 
                }
             } else { # gradient failure
                warnstr<-paste("Gradient not computable after method ",method[i],sep='')
                if (ctrl$dowarn) warning(warnstr)
                if (ctrl$trace>0) cat(warnstr,"\n") 
             }
          } # end kkt test
          ans$xtimes <- time
          # Do we want more information saved?
          if (ctrl$trace>0) { 
		cat("Save results from method ",meth,"\n") 
	  	print(ans)
	  }
	  if (ctrl$trace>0) { cat("Assemble the answers\n") }
          if (ctrl$trace>0) { 
              print(ans.ret)
              cat("add solution:\n")
              print(c(ans$par, ans$value, ans$fevals, ans$gevals, ans$nitns,
                              ans$convcode, ans$kkt1, ans$kkt2, ans$xtimes) )
          }                  
          ans.ret[meth, ] <- c(ans$par, ans$value, ans$fevals, ans$gevals, ans$nitns,
                              ans$convcode, ans$kkt1, ans$kkt2, ans$xtimes)
          if (! gradOK) ngatend <- NA
          if (! hessOK) {
             nhatend <- NA
             hev <- NA
          }
      }  ## end post-processing of successful solution
      ans.details<-rbind(ans.details, list(method=meth, ngatend=ngatend, nhatend=nhatend, hev=hev, message=amsg))
      # 1303234 try list() not c()
      row.names(ans.details)[[i]]>=meth
      	if (ctrl$follow.on) {
		par <- ans$par # save parameters for next method
		if (i < nmeth && (ctrl$trace>0)) cat("FOLLOW ON!\n") # NOT trace ??
	}
    } ## end loop over method (index i)
    ansout <- NULL # default if no answers
    if (length(ans$par) > 0) { # cannot save if no answers
	ansout <- data.frame(ans.ret)# Don't seem to need drop=FALSE
        attr(ansout, "details")<-ans.details
    }
    ansout # return(ansout)
} ## end of optimx.run
