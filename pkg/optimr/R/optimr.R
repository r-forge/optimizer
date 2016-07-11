optimr <- function(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {

# Check if bounded
  bdmsk <- bmchk(par, lower=lower, upper=upper)
  control$have.bounds <- bdmsk$bounds # and set a control value
#  cat("control$have.bounds =",control$have.bounds,"\n")


  orig.method <- method
  orig.gr <- gr
  orig.fn <- fn
  savehess <- hessian
  hessian <- FALSE

  if (is.null(method)) method <- "Nelder-Mead"

  if (is.null(control$trace)) control$trace <- 0
  npar <- length(par)
  defctrl <- ctrldefault(npar) # could leave this out in most cases

  ## ?? check if maximize and fnscale conflict ??

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

  if (is.null(gr)) gr <- defctrl$defgrapprox
  if (is.character(gr)) {
     egr <- function(spar, ...){
        if (control$trace > 1) {
           cat("fnscale =",fnscale,"  pscale=")
           print(pscale)
           cat("gr:")
           print(gr)
           par <- spar*pscale
           cat("par:")
           print(par)
        }
        result <- do.call(gr, list(par, userfn=fn, ...)) * fnscale
     }
  } else { 
    egr <- function(spar, ...) {
       par <- spar*pscale
       result <- gr(par, ...) * pscale * fnscale
    }
  } # end egr definition

  nlmfn <- function(spar, ...){
     f <- efn(spar, ...)
     g <- egr(spar, ...)
     attr(f,"gradient") <- g
     attr(f,"hessian") <- NULL # ?? maybe change later
     f
  }
# ?? do we want ehess ?    Not at 150714

## Masks 
   maskmeth <- control$maskmeth
   bdmsk <- bdmsk$bdmsk # Only need the masks bit from here on
   if (any(bdmsk == 0) ) {
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
      if (method== "Nelder-Mead" || 
          method == "BFGS" || 
          method == "L-BFGS-B" || 
          method == "CG" || 
          method == "SANN") {
        # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
        mcontrol$maxit <- control$maxit 
        mcontrol$trace <- control$trace
	mcontrol$parscale <- control$parscale # Use internal scaling
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
              ans <- try(optim(par=par, fn=orig.fn, gr=orig.gr, 
                      lower=lower, upper=upper, method="L-BFGS-B", hessian=FALSE, 
                       control=mcontrol, ...))
          }
        } else {
#          cat("calling optim() with no bounds\n")
          ans <- try(optim(par=par, fn=orig.fn, gr=orig.gr, 
                method=method, hessian=FALSE, control=mcontrol, ...))
#          print(ans)

        }
        if (class(ans)[1] == "try-error") { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
                errmsg <- "optim method failure\n"
                if (method != "L-BFGS-B") errmsg <- paste("optim() with bounds ONLY uses L-BFGS-B: ", errmsg)
		if (control$trace>0) cat(errmsg)
		ans$value <- control$badval
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
		ans$value <- control$badval
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
        if (is.null(gr)) { stop("optimr -- nlm -- we do not allow gr = NULL") }
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
		ans$value <- control$badval
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
	if (control$have.bounds) { # 151220 -- this was not defined
   	   ans <- try(Rcgminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=bdmsk, control=mcontrol, ...))
	} else {
   	   ans <- try(Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
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
        ## return(ans)
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (method == "Rvmmin") { # Use Rvmmin routine (ignoring masks??)
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
	mcontrol$trace <- control$trace # 140902 Note no check on validity of values
	if (control$have.bounds) {
   	   ans <- try(Rvmminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=bdmsk, control=mcontrol, ...))
	} else {
   	   ans <- try(Rvmminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
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
        ## return(ans)
      }  ## end if using Rvmmin
## --------------------------------------------
      else 
      if (method == "hjn") {# Use JN Hooke and Jeeves
        if (control$trace > 1) cat("hjn\n")
        if (control$trace > 0) {
           cat("control$have.bounds =",control$have.bounds,"\n")
        }
        ans <- try(hjn(spar, efn, lower=slower, upper=supper, bdmsk=bdmsk, 
                        control=control, ...))
        if (class(ans)[1] != "try-error") {
            ## Need to check these carefully??
            ans$par <- ans$par*pscale
            ans$value <- ans$value*fnscale
            ans$message <- NA # Should add a msg ??
         } else {
            if (control$trace>0) cat("hjn failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         ## return(ans)
      }  ## end if using lbfgs
## --------------------------------------------
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD:", method, sep='')
             stop(errmsg, call.=FALSE)
      }
# Exit from routine
      ## optexit -- function for return from routine adding in hessian
      if (savehess) { # compute hessian
         if (is.null(orig.gr)) {
            hess <- hessian(orig.fn, ans$par, ...) # from numDeriv
         } else { 
            hess <- jacobian(orig.gr, ans$par, ...) # use Jacobian of gradient
         }
         ans$hessian <- hess
      }
      ans # last statement of routine
} ## end of optimr
