optimr <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {

## 180706: Problem with maxit. Likely issue is that opm gets ctrldefault and
## updates with actual control list. But here in nlm, maxit defaults to 100.
## 181130: same with lbfgsb3 ??
  npar <- length(par)
  ctrl <- ctrldefault(npar)
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
  }
  ## 180706: Should we try to streamline?

  if (is.null(method)) method <- ctrl$defmethod

  outmethod <- checksolver(method, ctrl$allmeth, ctrl$allpkg) # there will only be one! 
  if (is.null(outmethod)) {
		if (ctrl$trace > 0) cat("Solver ",method," missing\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 8888 # failed in run
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- paste("Missing method ",method)
                ans$hessian <- NULL
                return(ans) # can't proceed without solver
  }

# Check if bounded
  bdmsk <- bmchk(par, lower=lower, upper=upper, shift2bound=TRUE)
  if (bdmsk$parchanged) warning("Parameter(s) changed to nearest bounds\n")
  ctrl$have.bounds <- bdmsk$bounds # and set a control value

  orig.method <- method
  orig.gr <- gr
  orig.fn <- fn

  if (is.null(hessian) ){
     savehess <- FALSE
  } else { savehess <- hessian } # logical -- whether to save hessian 

#  if (is.null(ctrl$trace)) ctrl$trace <- ctrl$trace

  if (is.null(ctrl$parscale)) { 
        pscale <- rep(1,npar)
        if(ctrl$trace > 0) { cat("Unit parameter scaling\n") }
  } else { 
        pscale <- ctrl$parscale 
        if(ctrl$trace > 0) {
          cat("Parameter scaling:")
          print(pscale)
        }
  }
  spar <- par/pscale # scaled parameters
  slower <- -Inf
  supper <- Inf # to ensure defined
  if (ctrl$have.bounds) {
    slower <- lower/pscale
    supper <- upper/pscale
  }
  fnscale <- 1 # default to ensure defined
  if (is.null(ctrl$fnscale)) {
     if (! is.null(ctrl$maximize) && ctrl$maximize ) {fnscale <- -1}
  else if (! is.null(ctrl$maximize)) {
          if ( (ctrl$fnscale < 0) && ctrl$maximize) {fnscale <- -1} # this is OK
          else stop("ctrl$fnscale and ctrl$maximize conflict")
       } # end ifelse
  } # end else
  ctrl$origmaximize <- ctrl$maximize # save the original in case we need it
  ctrl$maximize <- FALSE # and ensure we minimize
  ctrl$fnscale <- fnscale # to ensure set again

# 160615 -- decided to postpone adding nloptr

  efn <- function(spar, ...) {
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn(par, ...) * fnscale
  }

  appgr<-FALSE # so far assuming analytic gradient
#  DO NOT PROVIDE A DEFAULT -- LET METHOD DO THIS
#  if (is.null(gr)) gr <- ctrl$defgrapprox
#  if (is.null(gr)) cat("gr is NULL\n")
  if (is.character(gr)) {
     appgr <- TRUE # to inform us that we are using approximation
     egr <- function(spar, ...){
        if (ctrl$trace > 1) {
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

## cat("Now check if we need scaled hessian\n")
## cat("is.null(hess) = ",is.null(hess),"\n")

## if (! is.null(hess)) {
##    cat("Hessian at parameters:\n")
##    print(hess(par,...))
## } 

  if (is.null(hess)) { ehess <- NULL}
  else { ehess <- function(spar, ...) {
                      par <- spar*pscale
                      result <- hess(par, ...) * pscale * pscale * fnscale
                      result
                  }
  }

## cat("Is ehess present? is.null(ehess) =",is.null(ehess),"\n")
## if (! is.null(ehess)) {
##   cat("eHessian at scaled parameters:\n")
##   print(ehess(spar,...))
## }


  if (appgr && (ctrl$trace>0)) cat("Using numerical approximation '",gr,"' to gradient in optimru()\n")

  nlmfn <- function(spar, ...){
     f <- efn(spar, ...)
     if (is.null(egr)) {g <- NULL} else {g <- egr(spar, ...)}
     attr(f,"gradient") <- g
     if (is.null(ehess)) { h <- NULL } else {h <- ehess(spar, ...)}
     attr(f,"hessian") <- h
     f
  }

## cat("Check nlmfn at spar\n")
## print(nlmfn(spar, ...))


## Masks 
   maskmeth <- ctrl$maskmeth
   msk <- bdmsk$bdmsk # Only need the masks bit from here on
   if (any(msk == 0) ) {
      if ( !(method %in% maskmeth) ) {
         stopmsg <- paste("Method ",method," cannot handle masked (fixed) parameters")
         stop(stopmsg)
      }
      if (ctrl$trace > 0) cat("Masks present\n")
   }

# replacement for optim to minimize using a single method

# time is in opm(), but not here
# The structure has   par, value, counts, convergence, message, hessian

# Run a single method

# expand bounds
  if (length(lower) == 1 && is.finite(lower) ) lower<-rep(lower,npar)
  if (length(upper) == 1 && is.finite(upper) ) upper<-rep(upper,npar)

  mcontrol <- list() # define the control list

# Methods from optim()
  if (method == "Nelder-Mead" || 
      method == "BFGS" || 
      method == "L-BFGS-B" || 
      method == "CG" || 
      method == "SANN") {
      # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
      mcontrol$maxit <- ctrl$maxit 
      if (! is.null(ctrl$maxit)) {mcontrol$maxit <- ctrl$maxit}
      mcontrol$trace <- ctrl$trace
      mcontrol$parscale <- NULL # using user fn 
      mcontrol$fnscale <- NULL
##      mcontrol$fnscale <- ctrl$fnscale # 180313 Carlo Lapid ?? wrong, use efn, egr

# Note: hessian always FALSE in these calls. But savehess may recover it.

#        cat("Before optim() call - ctrl$have.bounds =",ctrl$have.bounds,"\n")
      if (ctrl$have.bounds) {
        if (method != "L-BFGS-B") {
            errmsg <- "optim() can only handle bounds with L-BFGS-B\n"
            if (ctrl$trace > 0) cat(errmsg,"\n")
            ans <- list()
            class(ans)[1] <- "try-error"
            warning("optimr: optim() with bounds ONLY uses L-BFGS-B")
        } else {
            ans <- try(optim(par=par, fn=efn, gr=egr, 
                      lower=lower, upper=upper, method="L-BFGS-B", hessian=FALSE, 
                       control=mcontrol, ...))
          }
        } else {
          ans <- try(optim(par=par, fn=efn, gr=egr, 
                method=method, hessian=FALSE, control=mcontrol, ...))
        }
        if (class(ans)[1] == "try-error") { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
                errmsg <- "optim method failure\n"
                if (method != "L-BFGS-B") 
                     errmsg <- paste("optim() with bounds ONLY uses L-BFGS-B: ", errmsg)
		if (ctrl$trace>0) cat(errmsg)
		ans$value <- ctrl$badval
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
        mcontrol$eval.max <- ctrl$maxfeval
	if ( is.null(ctrl$trace) || is.na(ctrl$trace) || ctrl$trace == 0) { 
		mcontrol$trace = 0
	} else { 
		mcontrol$trace = 1 # this is EVERY iteration. nlminb trace is freq of reporting.
	}
        ans <- try(nlminb(start=spar, objective=efn, gradient=egr, hessian=ehess, lower=slower, 
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
		if (ctrl$trace>0) cat("nlminb failure\n")
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- "nlminb failed" # 180318 change from NULL
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using nlminb
## --------------------------------------------
      else if (method == "nlm") { # Use stats package nlm routine
#        if (is.null(gr)) { stop("optimr -- nlm -- we do not allow gr = NULL") }
	if (! is.null(ctrl$maxit) ) {iterlim <- ctrl$maxit }
        else { iterlim <- 100 }
	print.level <- 0 
        errmsg <- NULL
        if (ctrl$have.bounds) {
              if(ctrl$trace > 0) cat("nlm cannot handle bounds\n")
              errmsg <- "nlm cannot handle bounds\n"
            ##  stop("nlm tried with bounds")
            ans <- list()
            class(ans)[1] <- "try-error"
        } else {
          if (! is.null(ctrl$trace) && (ctrl$trace > 0) ) {print.level <- 2 } 
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
                ans$message <- paste("nlm: Convergence indicator (code) = ",ans$code)
                ans$code <- NULL
	} else {
		if (ctrl$trace > 0) cat("nlm failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- "nlm failed" # 180318 change from NULL
                ans$hessian <- NULL
        }
        print.level <- NULL # clean up
        ## return(ans)
      } # end if using nlm
## --------------------------------------------
      else if (method == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
        mcontrol$trace <- ctrl$trace
        mcontrol$maxit <- ctrl$maxit # 151217 JN
        if (! is.null(egr)) {
  	  if (ctrl$have.bounds) { # 151220 -- this was not defined
            # 170919 -- explicit reference to package
   	    ans <- try(Rcgminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
	  } else {
   	     ans <- try(Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
	  }
        }
        if (!is.null(egr) && (class(ans)[1] != "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (ctrl$trace>0) {
                    cat("Rcgmin failed for current problem \n")
                    if(is.null(egr)) cat("Note: Rcgmin needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
      else if (method == "Rcgmin2") { # Use Rcgmin2 routine (ignoring masks)
        mcontrol$trace <- ctrl$trace
        mcontrol$maxit <- ctrl$maxit # 151217 JN
        if (! is.null(egr)) {
  	  if (ctrl$have.bounds) { # 151220 -- this was not defined
            # 170919 -- explicit reference to package
   	    ans <- try(Rcgmin2::Rcgminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
	  } else {
   	     ans <- try(Rcgmin2::Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
	  }
        }
        if (!is.null(egr) && (class(ans)[1] != "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (ctrl$trace>0) {
                    cat("Rcgmin2 failed for current problem \n")
                    if(is.null(egr)) cat("Note: Rcgmin2 needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                if(is.null(egr)) {
                   ans$message <- "Must specify gradient function for Rcgmin2"
                   ans$convergence <- 9998 # for no gradient where needed
                }
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rcgmin2
## --------------------------------------------
      else if (method == "Rvmmin") { # Use Rvmmin routine (ignoring masks??)
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval
	mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values
	if (! is.null(egr)) {
          ans <- try(Rvmmin(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
        }
        if (ctrl$trace > 2) {
            cat("Rvmmin ans:")
            print(ans)
        }
        if (! is.null(egr) && (class(ans)[1] != "try-error")) {
            ans$par <- ans$par*pscale
            ans$bdmsk <- NULL
        } else {
            if (ctrl$trace>0) {
                cat("Rvmmin failed for current problem \n")
                if(is.null(egr)) cat("Note: Rvmmin needs gradient function specified\n")
            }
	    ans<-list() # ans not yet defined, so set as list

	    ans$value <- ctrl$badval
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
      else if (method == "Rvmminq") { # Use Rvmminq routine (ignoring masks??)
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval
	mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values
	if (! is.null(egr)) {
          ans <- try(Rvmminq(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
        }
        if (ctrl$trace > 2) {
            cat("Rvmminq ans:")
            print(ans)
        }
        if (! is.null(egr) && (class(ans)[1] != "try-error")) {
            ans$par <- ans$par*pscale
            ans$bdmsk <- NULL
        } else {
            if (ctrl$trace>0) {
                cat("Rvmminq failed for current problem \n")
                if(is.null(egr)) cat("Note: Rvmminq needs gradient function specified\n")
            }
	    ans<-list() # ans not yet defined, so set as list

	    ans$value <- ctrl$badval
	    ans$par<-rep(NA,npar)
	    ans$counts[1] <- NA # save function and gradient count information
	    ans$counts[2] <- NA # save function and gradient count information
	    ans$message <- NULL        
            if(is.null(egr)) {
               ans$message <- "Must specify gradient function for Rvmminq"
               ans$convergence <- 9998 # for no gradient where needed
            }
            ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rvmminq
## --------------------------------------------
      else if (method == "snewton") { # Use snewton routine (no bounds or masks??)
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval # changed from maxfevals 180321
       	mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values
       	ans<-list(par=NA, value=NA, counts=NA, convergence=-1,
                  message=NA, hessian=NULL) # ans not yet defined, so set as list
       	if ( (! is.null(egr)) && (! is.null(ehess)) ) {
          if (ctrl$have.bounds) { # 170919 make package explicit
           stop("snewton does not handle bounds") # 
          } 
       	  tans <- try( snewton(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol,...))
          if (ctrl$trace>1) {
              cat("snewton returns tans:")
              print(tans)
          }
          if  (class(tans)[1] == "try-error") { 
             ans$message <- "snewton failed"
             ans$convergence <- 9999
             if (ctrl$trace>0) {
               cat(ans$message,"\n")
             }
          }
       	} else {
       	   if(is.null(egr)) {
       	     ans$message <- "Must specify gradient function for snewton"
       	     ans$convergence <- 9998 # for no gradient where needed
       	     warning("Note: snewton needs gradient function specified")
       	   }
       	   if(is.null(ehess)) {
       	     ans$message <- "Must specify Hessian function (hess) for snewton"
       	     ans$convergence <- 9997 # for no gradient where needed
       	     warning("Note: snewton needs Hessian function (hess) specified")
       	   }
       	}
       	if (ans$convergence > 9996){
       	       ans$value <- ctrl$badval
               ans$par<-rep(NA,npar)
               ans$counts[1] <- NA # save function and gradient count information
               ans$counts[2] <- NA # save function and gradient count information
               # Note: in optim() no provision for hessian count
               ans$hessian <- NULL
               if (ctrl$trace>1) {
                  cat("snewton falure ans:")
                  print(ans)
               }
        } else { # have an answer
              ans$par <- tans$par*pscale
              ans$value <- tans$value
              attr(ans, "gradient") <- tans$grad
              if(hessian) ans$hessian <- tans$Hess
              ans$counts[1] <- tans$counts$nfn
              ans$counts[2] <-  tans$counts$ngr
              ans$message <- tans$message 
              ans$convergence <- tans$convcode
              tans <- NULL # probably unnecessary, but for safety
              if (ctrl$trace>1) {
                   cat("rejigged ans:")
                   print(ans)
              }
             } # end have answer
         ## return(ans)
      }  ## end if using snewton
  ## --------------------------------------------
  else if (method == "snewtonm") { # Use snewtonm routine (no bounds or masks??)
    mcontrol$maxit <- ctrl$maxit
    mcontrol$maxfeval <- ctrl$maxfeval # changed from maxfevals 180321
    mcontrol$trace <- ctrl$trace # 140902 Note no check on validity of values
    ans<-list(par=NA, value=NA, counts=NA, convergence=-1,
                  message=NA, hessian=NULL) # ans not yet defined, so set as list
    if ( (! is.null(egr)) && (! is.null(ehess)) ) {
      if (ctrl$have.bounds) { # 170919 make package explicit
        stop("snewtonm does not handle bounds") # ?? error -- doesn't handle bounds
      } 
      tans <- try( snewtonm(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol,...))
      if (ctrl$trace>0) {
           cat("snewtonm returns tans:")
           print(tans)
      }
      if  (class(tans)[1] == "try-error") { 
        ans$message <- "snewtonm failed"
        ans$convergence <- 9999
        if (ctrl$trace>0) {
          cat(ans$message,"\n")
        }
      }
    } else {
      if(is.null(egr)) {
        ans$message <- "Must specify gradient function for snewtonm"
        ans$convergence <- 9998 # for no gradient where needed
        warning("Note: snewtonm needs gradient function specified")
      }
      if(is.null(ehess)) {
        ans$message <- "Must specify Hessian function (hess) for snewtonm"
        ans$convergence <- 9997 # for no Hessian where needed
        warning("Note: snewtonm needs Hessian function specified")
      }
    } # end of fails
    if (ans$convergence > 9996){ # Bad solution
       ans$value <- ctrl$badval
       ans$par<-rep(NA,npar)
       ans$counts[1] <- NA # save function and gradient count information
       ans$counts[2] <- NA # save function and gradient count information
       # Note: in optim() no provision for hessian count
       ans$hessian <- NULL
       if (ctrl$trace>1) {
          cat("snewton falure ans:")
          print(ans)
       }
    } else { # have an answer
#      cat("copy answer with tans$par:\n") 
#      print(tans$par)              
      ans$par <- tans$par*pscale
      ans$value <- tans$value
      attr(ans, "gradient") <- tans$grad
      if(hessian) ans$hessian <- tans$Hess
      ans$counts[1] <- tans$counts$nfn
      ans$counts[2] <-  tans$counts$ngr
      ans$message <- tans$message 
      ans$convergence <- tans$convcode
      tans <- NULL # probably unnecessary, but for safety
      if (ctrl$trace>1) {
         cat("rejigged ans:")
         print(ans)
      }
    } # end have answer
    ans
  }  ## end if using snewtonm
  ## --------------------------------------------
      else if (method == "hjn") {# Use JN Hooke and Jeeves
        if (ctrl$trace > 0) { 
           # this function is in optimr, so does not need explicit package
           cat("hjn:ctrl$have.bounds =",ctrl$have.bounds,"\n")
           cat("optimr - hjn - msk:")
           print(msk)
        }
        # 180327 Cannot maximize with hjn itself.
        mcontrol <- ctrl # copy
        mcontrol$maximize <- NULL # and null out maximize
        ans <- try(hjn(spar, efn, lower=slower, upper=supper, bdmsk=msk, 
                        control=mcontrol, ...))
        if (class(ans)[1] != "try-error") {
            ## Need to check these carefully??
            ans$par <- ans$par*pscale
            ans$value <- ans$value*fnscale
            ans$message <- NA # Should add a msg ??
         } else {
            if (ctrl$trace > 0) cat("hjn failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- ctrl$badval
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
        mcontrol$maxit <- ctrl$maxit
        mcontrol$maxfeval <- ctrl$maxfeval
        if (ctrl$trace > 0) { 
            mcontrol$trace <- TRUE
            if (ctrl$trace > 1) mcontrol$triter <- 1 # default is 10
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
		if (ctrl$trace > 0) cat("spg failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
        if (is.null(ctrl$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- ctrl$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
        errmsg <- NULL
        if (ctrl$have.bounds) {
              if (ctrl$trace > 0) cat("ucminf cannot handle bounds\n")
              errmsg <- "ucminf cannot handle bounds\n"
              stop(errmsg)
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
		        if (ans$convergence == 1 
		            || ans$convergence == 2 
		            || ans$convergence == 4) {
         		        ans$convergence <- 0
		        } 
            ans$par <- ans$par*pscale
        	  ans$counts[1] <- ans$info[4]
        	  ans$counts[2] <- ans$info[4] # calls fn and gr together
        	  ans$info <- NULL # to erase conflicting name
        	  ans$nitns <- NULL
            ans$hessian <- NULL
            ans$invhessian.lt <- NULL
		        if (ctrl$trace > 0) cat("ucminf message:",ans$message,"\n")
            } else { # ucminf failed
            		if (ctrl$trace > 0) cat("ucminf failed for this problem\n")
		            ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		            ans$value <- ctrl$badval
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
	if (ctrl$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}
	ans<-list() # ans not yet defined, so set as list
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
        if (is.null(egr)) { ## fixed msg below (referred to lbfgs) 170214
            if (ctrl$trace > 0) cat("Rtnmin MUST have gradient provided\n")
            errmsg <- "Rtnmin MUST have gradient provided"
            class(ans)[1] <- "try-error"            
        } else {
           if (ctrl$have.bounds) {
   	      ans <- try(tnbc(x=spar, fgfun=nlmfn, lower=slower,
                   upper=supper, trace=mcontrol$trace, ...))
           } else {
   	      ans <- try(tn(x=spar, fgfun=nlmfn, trace=mcontrol$trace, ...))
	   }
        }
        if (class(ans)[1] == "try-error") {
        	if (ctrl$trace>0) cat("Rtnmin failed for current problem \n")
                ans$convergence <- 9999 # failed in run
	        ans$message <- "Rtnmin failed fo current problem"        
                if (is.null(egr)) {
                   ans$convergence <- 9998
                   ans$message <- errmsg
                   ans$value <- 1234567E20
                } 
		ans$value <- ctrl$badval
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
  	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace
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
		if (ctrl$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        if (ctrl$have.bounds) {
            warning("Cannot use uobyqa with bounds")
		if (ctrl$trace > 0) cat("Cannot use uobyqa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
		if (ctrl$trace > 0) cat("uobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
        if (ctrl$trace > 1) cat("Trying newuoa\n")
	mcontrol$maxfun <- ctrl$maxfeval
        mcontrol$iprint <- ctrl$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        if (ctrl$have.bounds) {
            warning("Cannot use newuoa with bounds")
		if (ctrl$trace > 0) cat("Cannot use newuoa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
		if (ctrl$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- ctrl$badval
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
           if (ctrl$trace>0) cat("nmkb cannot start if on any bound \n")
           warning("nmkb() cannot be started if any parameter on a bound")
           ans <- list() # ans not yet defined, so set as list
           ans$value <- ctrl$badval
           ans$par <- rep(NA,npar)
           ans$convergence <- 9999 # failed in run - ?? consider special code for nmkb on bounds
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
           ans$hessian <- NULL
        } else { # ok to proceed with nmkb()
        if (! is.null(ctrl$maxit)) { 
	   mcontrol$maxfeval <- ctrl$maxit
	} else {
	   mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?
	}
         if (ctrl$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         if (ctrl$have.bounds) {
            ans <- try(dfoptim::nmkb(par=spar, fn=efn, lower = slower, 
              upper = supper, control=mcontrol, ...))
         } else {# 170919 explicit package in call
            ans <- try(dfoptim::nmk(par=spar, fn=efn, control=mcontrol, ...))
         }
         if (ctrl$trace > 1) {
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
           if (ctrl$trace>0) cat("nmkb failed for current problem \n")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- ctrl$badval
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
         if (ctrl$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
         mcontrol$maxfeval <- ctrl$maxfeval
         if (ctrl$have.bounds) {
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
            if (ctrl$trace>0) cat("hjkb failed for current problem \n")
            ans <- list(value=ctrl$badval, par=rep(NA,npar), message="Failed",
                convergence=9999)
            ans$counts[1]<- NA
            ans$counts[2]<- NA 
            ans$hessian <- NULL
         }
         ## return(ans)
      }  ## end if using hjkb
## --------------------------------------------
      else if (method == "lbfgsb3") {# Use 2011 L-BFGS-B wrapper
        if (ctrl$trace > 1) cat("lbfgsb3\n")
        mcontrol$trace <- ctrl$trace
        mcontrol$maxit <- ctrl$maxit
        if (ctrl$trace > 2) cat("lbfgsb3: maxit=",mcontrol$maxit,"\n")
# 170924 no longer needed
##        if (ctrl$trace < 1) {mcontrol$iprint <- -1} else {mcontrol$iprint <- ctrl$trace} 
        if (ctrl$trace > 2) cat("lbfgsb3:control$have.bounds =",ctrl$have.bounds,"\n")
        if (ctrl$have.bounds) { ## Note call uses prm not par
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
            if (ctrl$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- ctrl$badval
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
      else if (method == "lbfgsb3c") {# Use 2011 L-BFGS-B wrapper
        if (ctrl$trace > 1) cat("lbfgsb3c\n")
        mcontrol$maxit <- ctrl$maxit
        if (ctrl$trace > 3) cat("lbfgsb3c: maxit=",mcontrol$maxit,"\n")
        mcontrol$trace <- ctrl$trace
# 170924 no longer needed
##        if (ctrl$trace < 1) {mcontrol$iprint <- -1} else {mcontrol$iprint <- ctrl$trace} 
        if (ctrl$trace > 0) cat("lbfgsb3c:control$have.bounds =",ctrl$have.bounds,"\n")
        if (ctrl$have.bounds) { ## Note call uses prm not par
            slower <- lower/pscale
            supper <- upper/pscale
            ans <- try(lbfgsb3c::lbfgsb3c(par=spar, fn=efn, gr=egr, lower = slower, 
                upper = supper, control=mcontrol, ...)) # explicit pkg in call 170919
        } else {
            ans <- try(lbfgsb3c::lbfgsb3c(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
        }
        if (class(ans)[1] != "try-error") {
 ## Need to check these carefully??
#            ans$convergence <- 0
            ans$par <- ans$par*pscale
#            ans$prm <- NULL
#            ans$value<-as.numeric(ans$f)
#            ans$f <- NULL
#            ans$g <- NULL ## perhaps keep -- but how??
#            ans$hessian <- NULL
            ans$message <- NA
            ans$niter <- NULL # loss of information
         } else {
            if (ctrl$trace>0) cat("lbfgsb3 failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- ctrl$badval
            ans$par<-rep(NA,npar)
            ans$convergence<-9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         ## return(ans)
      }  ## end if using lbfgsb3c
## --------------------------------------------
      else if (method == "lbfgs") {# Use unconstrained method from lbfgs package
        if (ctrl$trace > 1) cat("lbfgs\n")
        if (ctrl$trace < 1) {invisible <- 1} else {invisible <- 0}
        if (ctrl$trace > 1) cat("lbfgs:control$have.bounds =",ctrl$have.bounds,"\n")
        ans <- list() # to define the answer object
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
##      cat("in lbfgs section, ctrl$have.bounds=",ctrl$have.bounds,"\n")
        if (ctrl$have.bounds) {
              cat("ctrl$have.bounds seems TRUE\n")
              if (ctrl$trace > 0) cat("lbfgs::lbfgs cannot handle bounds\n")
              errmsg <- "lbfgs::lbfgs cannot handle bounds\n"
            ##  stop("lbfgs::lbfgs tried with bounds")
            class(ans)[1] <- "try-error"            
        }
        if (is.null(egr)) {
            if (ctrl$trace > 0) cat("lbfgs::lbfgs MUST have gradient provided\n")
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
            if (ctrl$trace>0) cat("lbfgs failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- ctrl$badval
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
    if (ctrl$trace > 1) cat("subplex\n")
    if (ctrl$trace < 1) {invisible <- 1} else {invisible <- 0}
    if (ctrl$trace > 1) cat("subplex:control$have.bounds =",ctrl$have.bounds,"\n")
    ans <- list() # to define the answer object
    if (ctrl$trace > 0) warning("subplex has no trace mechanism")
    class(ans)[1] <- "undefined" # initial setting
    if (ctrl$have.bounds) {
      cat("ctrl$have.bounds seems TRUE\n")
      if (ctrl$trace > 0) cat("subplex::subplex cannot handle bounds\n")
      stop("subplex::subplex cannot handle bounds")
    }
    if (class(ans)[1] == "undefined"){
       ans <- try(subplex::subplex(par=spar, fn=efn, control=list(maxit=ctrl$maxfeval)))
    }
    if ((class(ans)[1] != "try-error") && (ans$convergence != -2)) {
       ## Need to check these carefully??
       ans$par <- ans$par*pscale
       ans$value <- ans$value*fnscale
       ans$counts[1] <- ans$count
       ans$counts[2] <- NA
       ans$count <- NULL
       ccode <- ans$convergence
       ans$convergence <- 9999
       ans$message <- paste("subplex:",ans$message)
       if ((ccode == 0) || (ccode == 1)) {
                ans$convergence <- 0 
                if (ccode == 0) { ans$message <- "subplex: success" }
                else { ans$message <- "subplex: Limit of precision reached" }
           } # converged OK
           else {if (ccode == -1) {
                    ans$convergence <- 1
                    ans$message <- "subplex: function evaluation limit reached"
                } # effort limit
       }
    } else { 
      if (ccode == -2) {
            ans$convergence <- 20
      }
      else { ans$convergence <- 9999}
      ans$value <- ctrl$badval
      ans$par <- rep(NA,npar)
      ans$counts[1] <- NA
      ans$counts[2] <- NA
      ans$hessian <- NULL
    }
  }  ## end if using subplex
## --------------------------------------------
## END OF optimrx extra methods
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD:", method, sep='')
             stop(errmsg, call.=FALSE)
      }
# Exit from routine
      ans$value <- ans$value * fnscale # reset for maximum
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
