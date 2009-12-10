optimx <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=NULL,
            control=list(),
             ...) {
##### OPEN ISSUES: (any date order)
# 091018 -- ?? new control "starttest" so we can skip them
# 091018 -- ?? masks?
# 091018 -- ?? exit when initial point is infeasible -- now just exit, should move to next problem?
# 090929 -- ?? control structure is set up badly -- needs simplification 
# 090729 -- ?? simplify structure of answers -- avoid repetition and make easier to access
# 090531 -- ?? need SNewton back in to keep Ramsay and Co. on side
# 090601 -- ?? Do we want hessevals to count Hessian evaluations?
# 090612 -- ?? There may be better choices for the tolerances in tests of equality for gradient / kkt tests.

##### IMPLEMENTED: (reverse date order)
# 091026 -- Remove SANN because no real conv indicator
# 091018 - 090531 Use function code joint with funtest and funcheck in initial checks -- not fully possible
# 091018 -- decided not to Put optimx.dev into package to provide for local source packages etc.
# 090923 -- add bobyqa (from minqa)
# 090923 -- add Rcgmin for large scale problems
# 090729 -- put hessian argument back as an alternate for kkt
# 090729 -- should have kkt only carried out for n<50 (no gradient), n<500 (with gradient and Jacobian
#            of gradient to get Hessian)
# 090531 -- KKT stuff
# 090531 -- control for keeping failed runs (save.failures)
# 090531 Decided to omit 'trust' from package methods for the time being.
# 090527 -- added follow.on
# 090511 What should be default method(s)? 
#         090601: to provide compatibility with optim(), Nelder-Mead is put first. 
#	          The choice of BFGS second is open to discussion. JN


#  A wrapper function to integrate major optimization packages in R
#  For now, we only consider algorithms for unconstrained and box-constrained optimization
#  This function can implement "multiple" algorithms
#
# Input:
#  par = a single vector of starting values
#  fn = objective function (assumed to be sufficeintly differentiable)
#  gr = name of a function to compute the (analytic) gradient of the objective function
#  hess = name of a function to compute the (analytic) Hessian of the objective function
#         Editorial note: Capitalize Hessian after the name of Otto Hesse. JN
#  method = a character vector denoting all the algorithms to be executed (in the specified order)
#      Complete names are not needed. A partial matching is attempted.
#  hessian = logical variable that, if present, is equivalent to control$kkt. If TRUE, it causes
#      optimx to try to compute an approximation to the Hessian matrix at the final set of parameters.
#  control = list of control information, in particular
#      trace = an integer controlling output (note different methods may use logicals
#         trace = 0 gives no output, positive numbers give more output for greater values
#      follow.on = TRUE or FALSE. If TRUE, and there are multiple methods, then the last set of 
#         parameters from one method is used as the starting set for the next. 
#      save.failures = TRUE if we wish to keep "answers" from runs where the method does not 
#         return conv==0. FALSE otherwise (default).
#      maximize = TRUE if we want to maximize rather than minimize a function. (Default FALSE)
#         090601: Not yet implemented for nlm, nlminb, ucminf. However, there is a check to avoid
#                 usage of these codes when maximize is TRUE.
#      all.methods = TRUE if we want to use all available (and suitable) methods
#      starttests = TRUE assumes we do initial tests for feasibility
#
#
#
# Output:
# ans = an object containing two sets of information:
# essential output = a data frame of converged parameters, function value at convergence, 
#    name of algorithm, convergence indicator, and function evaluation and gradient evaluation counts
# detailed output = this is a list with one component for each algorithm, and contains parameters, 
#    function values, convergence code, number of function and gradient evals, numerical gradient 
#    and hessian at convergence, eigenvalues of that hessian approximation, cpu time, and other 
#    information returned by algorithms, and name of the algorithm.
# detailed output can be accessed via the attribute called `details'
#
#  Authors:  Ravi Varadhan & John Nash
#  Date:  February 17, 2008
#  Changes: Ravi Varadhan - Date: May 29, 2009, John Nash - Latest: Oct 18, 2009
#
#################################################################
# Get real name of function to be minimized
  fname<-deparse(substitute(fn))
  if (!is.null(control$trace) && control$trace>0) {
	cat("fn is ",fname,"\n")
	tmp<-readline("Continue?")
  }

## Code more or less common to funtest, funcheck and optimx <<<
# Check parameters are in right form
  if(!is.null(dim(par))) stop("Parameter should be a vector, not a matrix!", call. = FALSE)
  if (! is.vector(par) ) {
	stop("The parameters are NOT in a vector")
  }
  npar<-length(par)
# Set control defaults
    ctrl <- list(
       follow.on=FALSE, 
       save.failures=FALSE, 
       trace=0, 
       maximize=FALSE,
       sort.result=TRUE,
       kkt=TRUE,
       all.methods=FALSE,
       starttests=TRUE
    ) # for now turn off sorting until we can work out how to do it nicely
    
# Note that we do NOT want to check on the names, because we may introduce 
#    new names in the control lists of added methods
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
# However, we do want to substitute the appropriate information. 
# hessian control gets copied to kkt control
    if (!is.null(hessian)){
	control$kkt<-hessian # Note: in optim, hessian is NOT in control settings.
    }
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
       if (onename %in% nctrl) {
           ctrl[onename]<-control[onename]
       } else {
           ctrl[onename]<-control[onename]
       }
    }
    if (is.null(control$kkt)) { # default is to turn off kkt for large matrices
      if (is.null(gr)) { # no analytic gradient
         if (npar > 50) {
           ctrl$kkt=FALSE # too much work when large number of parameters
           if (ctrl$trace) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      } else {
         if (npar > 500) {
            ctrl$kkt=FALSE # too much work when large number of parameters, even with analytic gradient
            if (ctrl$trace) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      }
    } else { # kkt is set
      if (control$kkt) {
        if (is.null(gr)) {
           if (npar > 50) {
             if (ctrl$trace) warning("Computing hessian for gr NULL, npar > 50, can be slow\n")
           }
        } else {
           if (npar > 500) {
             if (ctrl$trace) warning("Computing hessian even with gradient code, npar > 500, can be slow\n")
           }
        }
      }
    }
# Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) { have.bounds<-TRUE # set this for convenience
  } else { have.bounds <- FALSE }

  if (ctrl$starttests) {
	# Check parameters in bounds (090601: As yet not dealing with masks ??)
	#  bdmsk<-as.vector(bdmset[k, ])

	infeasible<-FALSE
	if (ctrl$trace > 0) cat("Function has ",npar," arguments\n")
 
	if (have.bounds) {
    	  # Expand bounds to vectors if needed
    	  if (length(lower)==1 ) lower <- rep(lower, npar)
    	  if (length(upper)==1 ) upper <- rep(upper, npar)
    	  bstate<-vector(mode="character", length=npar)
    	  for (i in 1:npar) {
     	    if ( (lower[i]<=par[i]) && (par[i]<=upper[i])) {
    	      bstate[i]<-" In Bounds "
            } else { 
            #   if(bdmsk[i]!=0) {
                  infeasible<-TRUE
            #   }
                if (lower[i]>par[i]) {bstate[i]<-" Out of Bounds LOW" } else { bstate[i]<-" Out of Bounds HIGH " }
            } # end if in bounds
#           if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bdmsk[i],"   ",bstate,"\n")
            if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bstate,"\n")
          } # end of for loop over parameter vector elements
	  if(infeasible) { 
        	stop("Infeasible point, no further tests")
	  } 
  	} # end have.bounds
        # Check if function can be computed
        firsttry<-try(finit<-fn(par, ...) ) 
        # Note: This incurs one EXTRA function evaluation because optimx is a wrapper for other methods
        if (class(firsttry) == "try-error") {
    	   infeasible <- TRUE
           stop("Cannot evaluate function at initial parameters")
        }
        # Also check that it is returned as a scalar
#       cat("finit=",finit,"\n")
#       print(str(finit))
#       cat("is.vector(finit)=",is.vector(finit),"\n")
#       cat("is.list(finit)=",is.list(finit),"\n")
#       cat("is.array(finit)=",is.array(finit),"\n")
#       cat("is.matrix(finit)=",is.matrix(finit),"\n")
#       cat("is.numeric(finit)=",is.numeric(finit),"\n")
#       cat("length(finit)=",length(finit),"\n")
       if (!(is.vector(finit) && (length(finit)==1)) || is.list(finit) || 
            is.matrix(finit) || is.array(finit) || ! is.numeric(finit) ) {
           stop("Function provided is not returning a scalar number")
       }
       if (is.infinite(finit) || is.na(finit)) {
          stop("Function returned is infinite or NA (non-computable)")
       }
  }

  # Check that we have the functions we need
  ipkgs <- installed.packages()
  if ("numDeriv" %in% ipkgs[,1]) require(numDeriv, quietly=TRUE) 
  else stop("Install package `numDeriv'", call.=FALSE)


  if(ctrl$starttests) {
     if(! is.null(gr)){ # check gradient
       gname <- deparse(substitute(gr))
       if (ctrl$trace) cat("Analytic gradient from function ",gname,"\n\n")
          fval <- fn(par,...) 
          gn <- grad(func=fn, x=par,...) # 
          ga <- gr(par, ...)
          # Now test for equality (090612: ?? There may be better choices for the tolerances.
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(gn-ga))/(1 + abs(fval)) >= teps) stop("Gradient function might be wrong - check it! \n", call.=FALSE)
       } else if (ctrl$trace) cat("Analytic gradient not made available.\n")

       if(! is.null(hess)){ # check Hessian
          hname <- deparse(substitute(hess))
          if (ctrl$trace) cat("Analytic hessian from function ",hname,"\n\n")
          hn <- hessian(func=fn, x=par,...) # ?? should we use dotdat
          ha <- hess(par, ...)
          # Now test for equality
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
       } else if (ctrl$trace) cat("Analytic Hessian not made available.\n")
   }
## >>> End of code common to funtest, funcheck and optimx, with mods for local needs

  # Set up the vector to return answers
  ans.ret <- vector(mode="list")
  # mode= is not strictly required. length defaults to 0. This sets up our answer vector.

  # List of methods in base or stats, namely those in optim(), nlm(), nlminb()
#  bmeth <- c("BFGS", "CG", "Nelder-Mead", "SANN", "L-BFGS-B", "nlm", "nlminb", "DEoptim")
# Remove SANN because it has no termination for optimality, only a maxit count for
#    the maximum number of function evaluations; remove DEoptim for now -- not useful 
#    for smooth functions. Code left in for those who may need it.
  bmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb")
  # List of methods in packages. 
  pmeth <- c("spg", "ucminf", "Rcgmin", "bobyqa")
  allmeth <- c(bmeth, pmeth)
  nomaxmeth <- c("nlm", "nlminb", "ucminf") # methods not (yet) allowing maximize
  # Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) allmeth <- c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "bobyqa") 
  if (ctrl$maximize) {
     for (meth in allmeth) {
         if (meth %in% nomaxmeth) { 
             if(ctrl$trace>0) cat("Method ",meth," cannot maximize -- dropping\n")
             allmeth<-allmeth[! (meth == allmeth)]
         }
     }
  }
  if (ctrl$all.methods) {
	method<-allmeth
	if (ctrl$trace>0) {
		cat("all.methods is TRUE -- Using all available methods\n")
		print(method)
	}
  }

  # Partial matching of method string allowed
  method <- unique(match.arg(method, allmeth, several.ok=TRUE) )
  nmeth <- length(method) # number of methods requested

  ## Check that methods are indeed available and loaded
  for (i in 1:nmeth) {
     cmeth <- method[i]
     if (ctrl$trace > 0) cat("Looking for method = ",cmeth,"\n")
     if (! (cmeth %in% allmeth) ) {
         errmsg <- paste(cmeth," not found in any list of methods available")
         stop(errmsg, call.=FALSE)
     } # otherwise the method is available, and just needs to be loaded
  } # end check methods available
  if (ctrl$trace>1) {
    cat("Methods to be used:")
    print(method)
  }
# Now make sure methods loaded
  if(any(method == "spg")) {
	if ("BB" %in% ipkgs[,1]) require(BB, quietly=TRUE)
	else  stop("Package `BB' Not installed", call.=FALSE)
	}
  if(any(method == "ucminf")) { 
	if ("ucminf" %in% ipkgs[,1]) require(ucminf, quietly=TRUE)
	else  stop("Package `ucminf' Not installed", call.=FALSE)
	}
  if(any(method == "Rcgmin")) { 
	if ("Rcgmin" %in% ipkgs[,1]) require(Rcgmin, quietly=TRUE)
	else  stop("Package `Rcgmin' Not installed", call.=FALSE)
	}
  if(any(method == "bobyqa")) { 
	if ("minqa" %in% ipkgs[,1]) require(minqa, quietly=TRUE)
	else  stop("Package `minqa' (for bobyqa) Not installed", call.=FALSE)
	}
#  if(any(method == "DEoptim")) { # Code removed as DEoptim not part of current set of methods
#	if ("DEoptim" %in% ipkgs[,1]) require(DEoptim, quietly=TRUE)
#	else  stop("Package `DEoptim' Not installed", call.=FALSE)
#	}

# Run methods
  times <- rep(0, nmeth)  # figure out how many methods and reserve that many times to record.
  names(times) <- method  # And label these times appropriately with the method.
  j <- 0  # j will be the index to the answers, one per each method, starting parameter pair

  for (i in 1:nmeth) { # loop over the methods
      meth <- method[i] # extract the method name
      conv <- -1 # indicate that we have not yet converged
      if (ctrl$trace) cat("Method: ", meth, "\n") # display the method being used
      # Extract control information e.g., maxit
      if (! is.null(ctrl$maxit) ) {
        if (length(ctrl$maxit) == 1 ) {
          if(ctrl$trace>0) cat("setting meth.maxit to common value ",ctrl$maxit,"\n")
          meth.maxit <- ctrl$maxit # We use a single value for all tries
        } else {
          if (length(ctrl$maxit) != nmeth) stop("ctrl$maxit has wrong number of elements", call.=FALSE)
          if(ctrl$trace>0) cat("setting meth.maxit for choice ",i," = ",ctrl$maxit[i],"\n")
          meth.maxit <- ctrl$maxit[i] # for the case with separate methods
        }
      } else { meth.maxit<-NULL }
      # create local control list for a single method -- this is one of the key issues for optimx
      mcontrol<-ctrl
      mcontrol[which(names(mcontrol) == "maxit")] <- meth.maxit
      mcontrol$follow.on<-NULL # And make sure that controls NOT in method are deleted (nulled)
      mcontrol$save.failures<-NULL
      mcontrol$sort.result<-NULL
      mcontrol$kkt<-NULL
      mcontrol$starttests<-NULL
      mcontrol$all.methods<-NULL
# Handle maximization if necessary
      if (meth != 'spg') mcontrol$maximize<-NULL # only spg has explicit maximize
      if (ctrl$maximize) { # Try to sort out "maximize" control
         if (any(meth == nomaxmeth)) {
		stop("At the moment you must explicitly set up maximization in your function for ", meth)
         }
         if (! is.null(ctrl$fnscale)) { 
 		stop("Mixing maximize and fnscale is dangerous. Please correct.")
         } 
      }
# Methods from optim()
      if (meth=="Nelder-Mead" || meth == "BFGS" || meth == "L-BFGS-B" || meth == "CG" || meth == "SANN") {
#        if (meth == "SANN") mcontrol$maxit<-10000 # !! arbitrary for now
        # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, SANN and CG
        time <- system.time(ans <- try(optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, 
                method=meth, control=mcontrol, ...), silent=TRUE))[1]
        # The time is the index=1 element of the system.time for the process, 
        # which is a 'try()' of the regular optim() function
        if (class(ans) != "try-error") { 
		ans$conv <- ans$convergence
		#      convergence: An integer code. '0' indicates successful convergence.
	        ans$fevals<-ans$counts[1] # save function and gradient count information
	        ans$gevals<-ans$counts[2]
        	ans$counts<-NULL # and erase the counts element now data is saved
        	ans$niter<-NULL # not used, so delete
#                if ((meth=="SANN") && (ans$fevals>=mcontrol$maxit)) {
#			if(ctrl$trace) cat("SANN at limit!")
#			# SANN uses rather different conventions to other optim() tools.
#                        ans$conv<-1
#		}
	} else { # bad result -- what to do
		ans$fefals<-NA
                ans$conv<-9999 # failed in run
		ans$value<-NA
        }
      }   # end if using optim() methods
## --------------------------------------------
      else if (meth == "nlminb") {
        # Here we use portLib routine nlminb rather than optim as our minimizer
        mcontrol$iter.max<-mcontrol$maxit # different name for iteration limit in this routine
        mcontrol$maxit<-NULL
        time <- system.time(ans <- try(nlminb(start=par, objective=fn, gradient=gr, lower=lower, 
		upper=upper, control=mcontrol,  ...), silent=TRUE))[1]
        if (class(ans) != "try-error") ans$conv <- ans$convergence
        # Translate output to common format and names
        ans$value<-ans$objective
        ans$objective<-NULL
        ans$fevals<-ans$evaluations[1]
        ans$gevals<-ans$evaluations[2]
	ans$evaluations<-NULL # cleanup
        ans$niter<-ans$iterations
        ans$iterations<-NULL
      }  ## end if using nlminb
## --------------------------------------------
      else if (meth == "nlm") { # Use stats package nlm routine
        ufn <- fn # we don't want to change the user function object, so create a copy
        if (!is.null(gr)) {
	   attr(ufn, "gradient") <- gr(par, ...) 
           if (!is.null(hess)) attr(ufn, "hessian") <- hess(par, ...) # Only change attibute if gr defined too
        }
        time <- system.time(ans <- try(nlm(f=ufn, p=par, ..., print.level=ctrl$trace), silent=TRUE))[1]
        # print(ans)
        if (class(ans) != "try-error") {
		ans$conv <- ans$code
		if (ans$conv == 1 | ans$conv == 2 | ans$conv == 3) ans$conv <- 0
        	# Translate output to common format
  	        if (ans$conv == 0) {
			ans$value<-ans$minimum
                        ans$par<-ans$estimate
			ans$estimate<-NULL
			ans$minimum<-NULL
                }
        	ans$fevals<-NA
        	ans$gevals<-NA # ?? need to fix this somehow in nlm code
        	ans$niter<-ans$iterations
        	ans$iterations<-NULL
	} else {
		if (ctrl$trace > 0) cat("nlm failed for this problem\n")
		ans$value<-Inf # 091209 to fix sorting
        	ans$conv<-9999
        	ans$fevals<-NA
        	ans$gevals<-NA 
        	ans$niter<-NA
	}
      } # end if using nlm
## --------------------------------------------
      else if (meth == "spg") { # Use BB package routine spg as minimizer
        time <- system.time(ans <- try(spg(par=par, fn=fn, gr=gr, lower=lower, upper=upper,  
		control=list(maxit=1500, trace=FALSE), ...), silent=TRUE))[1]
        if (class(ans) != "try-error") { 
   	   ans$conv <- ans$convergence
           ans$fevals<-ans$feval
           ans$feval<-NULL # to erase conflicting name
           ans$gevals<-NA # ??fixup needed
           ans$niter<-ans$iter
           ans$iter<-NULL
        } else { # spg failed
		if (ctrl$trace > 0) cat("spg failed for this problem\n")
        	ans$conv<-9999
        	ans$fevals<-NA
        	ans$gevals<-NA 
        	ans$niter<-NA
        }
      }  # end if using spg
## --------------------------------------------
      else if (meth == "ucminf") {
        ## Use ucminf routine
        mcontrol$meval<-mcontrol$maxit # different name in this routine
        mcontrol$maxit<-NULL
#        time <- system.time(ans <- try(ucminf(par=par, fn=fn, gr=gr,  control=mcontrol, ...), silent=TRUE))[1]
        time <- system.time(ans <- try(ucminf(par=par, fn=fn, gr=gr,  ...), silent=TRUE))[1]
        if (class(ans) != "try-error") {
		ans$conv <- ans$convergence
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
		if (ans$conv == 1 | ans$conv == 2 | ans$conv == 4) {
         		ans$conv <- 0
		} else {
			ans$conv <- ans$convergence
		} # Termination criteria are tricky here!  How to determine successful convergence?
                ans$convergence<-NULL
        	ans$fevals<-ans$info[4]
        	ans$gevals<-ans$info[4] # calls fn and gr together
        	ans$info<-NULL # to erase conflicting name
        	ans$niter<-NULL
		if (ctrl$trace > 0) cat("ucminf message:",ans$message,"\n")
        } else { # spg failed
		if (ctrl$trace > 0) cat("spg failed for this problem\n")
        	ans$conv<-9999
        	ans$fevals<-NA
        	ans$gevals<-NA 
        	ans$niter<-NA
        }
      }  ## end if using ucminf
## --------------------------------------------
###### progress point #########
      else if (meth == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
	bdmsk<-rep(1,npar)
        time <- system.time(ans <- try(Rcgmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, bdmsk=bdmsk,...), silent=TRUE))[1]
        if (class(ans) != "try-error") {
		ans$conv <- ans$convergence
	        ans$fevals<-ans$counts[1]
	        ans$gevals<-ans$counts[2]
		ans$value<-ans$value 
		cat("Rcgmin fval=",ans$value,"\n")
        } else {
		if (ctrl$trace>0) cat("Rcgmin failed for current problem \n")
       		ans$conv<-9999
	        ans$fevals<-NA
	        ans$gevals<-NA
		ans$value<-NA
        }
        ans$niter<-NULL
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (meth == "bobyqa") {# Use bobyqa routine from minqa package
        if (! is.null(mcontrol$maxit)) { 
		mcontrol$maxfun<-mcontrol$maxit
	} else {
		mcontrol$maxfun<-500000 # ?? default at 091018, but should it be changed?!!
	}
        mcontrol$iprint<-mcontrol$trace
        time <- system.time(ans <- try(bobyqa(par=par, fn=fn, lower=lower, upper=upper, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans) != "try-error") {
		ans$conv <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$conv <- 1 # too many evaluations
                }
	        ans$fevals<-ans$feval
	        ans$gevals<-NA
		ans$value<-ans$fval 
        } else {
		if (ctrl$trace>0) cat("bobyqa failed for current problem \n")
       		ans$conv<-9999
	        ans$fevals<-NA
	        ans$gevals<-NA
		ans$value<-NA
        }
        ans$niter<-NULL
      }  ## end if using bobyqa
## --------------------------------------------
      else if (meth == "DEoptim") {# Use DEoptim package
        stop("This code for DEoptim has not been verified")
        if (! is.null(mcontrol$maxit)) { 
		mcontrol$maxfun<-mcontrol$maxit
	} else {
		mcontrol$maxfun<-500000 # ?? default at 091018, but should it be changed?!!
	}
        mcontrol$iprint<-mcontrol$trace
#        time <- system.time(ans <- try(DEoptim(fn=fn, lower=lower, upper=upper, control=mcontrol,...), silent=TRUE))[1]
        if (class(ans) != "try-error") {
		ans$conv <- 0
                if (ans$feval > mcontrol$maxfun) {
			ans$conv <- 1 # too many evaluations
                }
	        ans$fevals<-ans$feval
	        ans$gevals<-NA
		ans$value<-ans$fval 
        } else {
		if (ctrl$trace>0) cat("bobyqa failed for current problem \n")
       		ans$conv<-9999
	        ans$fevals<-NA
	        ans$gevals<-NA
		ans$value<-NA
        }
        ans$niter<-NULL
      }  ## end if using DEoptim
## --------------------------------------------
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD: ", meth, sep='')
             stop(errmsg, call.=FALSE)
           }
## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      times[i] <- times[i] + time # Accumulate time for a single method (in case called multiple times)
      if (ctrl$trace) { cat("Post processing for method ",meth,"\n") }
      if ( ctrl$save.failures || (ans$conv == 0) ){  # Save the solution if converged or directed to save
         j <- j + 1 ## increment the counter for (successful) method/start case
         if (ctrl$trace) cat("Successful convergence! \n")  ## inform user we've succeeded
         # Testing final solution. Use numDeriv to get gradient and Hessian; compute Hessian eigenvalues
         ans$kkt1<-NA
	 ans$kkt2<-NA
         # cat("Post processing -- ctrl$kkt = ",ctrl$kkt," ans$conv=",ans$conv,"\n")
         if(ctrl$kkt && (ans$conv != 9999)) { # need to avoid test when method failed
           ngatend <- if (is.null(gr)) grad(fn, ans$par) else gr(ans$par) # Gradient at solution
           nhatend <- if (is.null(gr)) hessian(fn, ans$par) else jacobian(gr,ans$par) # numerical hessian at "solution"
           # For bounds constraints, we need to "project" the gradient and Hessian
           bset<-sort(unique(c(which(ans$par<=lower), which(ans$par>=upper))))
           nbds<-length(bset) # number of elements nulled by bounds
  	   # Note: we assume that we are ON, not over boundary, but use <= and >=. No tolerance is used.
           ngatend[bset]<-0 # "project" the gradient
           nhatend[bset,] <-0
           nhatend[, bset] <-0 # and the Hessian
           ans$ngatend <- ngatend
           ans$nhatend <- nhatend
           ans$evnhatend <- eigen(nhatend)$values
	   # test results
	   teps <- (.Machine$double.eps)^(1/3) # What tolerance should be used?
           ans$kkt1<-(max(abs(ngatend)) <= teps*(1.0+abs(ans$value)) ) # ?? Is this sensible?
           # now look at Hessian
           negeig<-(ans$evnhatend[npar] <= (-1)*teps*(1.0+abs(ans$value)))
           evratio<-ans$evnhatend[npar-nbds]/ans$evnhatend[1]
	   # If non-positive definite, then there will be zero eigenvalues (from the projection)
	   # in the place of the "last" eigenvalue and we'll have singularity.
           # WARNING: Could have a weak minimum if semi-definite.
	   ans$kkt2<- (evratio > teps) && (! negeig) # May want to think some more about this.??
        } # end kkt test
        ans$systime <- time
        # Do we want more information saved?
        if (ctrl$trace) { cat("Save results from method ",meth,"\n") }
        ans.ret[[j]] <- ans  ## save the answer. [[]] indexes the CONTENTS of the list
        ans.ret[[j]]$method <- method[i] # and we tag on the method with the $ linker
      }  ## end post-processing of successful solution
      if (ctrl$trace) { cat("Check if follow.on from method ",meth,"\n") }
      if (ctrl$follow.on) par <- ans$par # save parameters for next method
    } ## end loop over method (index i)
    if (ctrl$trace) { cat("Assemble the answers\n") }
    attr(ans.ret, "CPU times (s)") <- times ## save the accumulated times 
    if (ctrl$trace) { cat("Done CPU times\n") }
    pars <- lapply(ans.ret, function(x) x$par)
    if (ctrl$trace) { cat("Done parameters\n") }
    vals <- lapply(ans.ret, function(x) x$value)
    if (ctrl$trace) { cat("Done value\n") }
    meths <- lapply(ans.ret, function(x) x$method)
    if (ctrl$trace) { cat("Done method\n") }
    fevals<- lapply(ans.ret, function(x) x$fevals)
    if (ctrl$trace) { cat("Done fevals\n") }
    gevals<- lapply(ans.ret, function(x) x$gevals)
    if (ctrl$trace) { cat("Done gevals\n") }
    nitns <-  lapply(ans.ret, function(x) x$niter)
    if (ctrl$trace) { cat("Done niter\n") }
    convcode<- lapply(ans.ret, function(x) x$conv)
    if (ctrl$trace) { cat("Done conv\n") }
    kkt1<- lapply(ans.ret, function(x) x$kkt1)
    if (ctrl$trace) { cat("Done kkt1\n") }
    kkt2<- lapply(ans.ret, function(x) x$kkt2)
    if (ctrl$trace) { cat("Done kkt2\n") }
    if (ctrl$trace) { cat("Consolidate ans\n") }
    ans <- data.frame(cbind(par=pars, fvalues=vals, method=meths, fns=fevals, grs=gevals, 
	itns=nitns, conv=convcode, KKT1=kkt1, KKT2=kkt2))
    if (ctrl$trace) { cat("Add details\n") }
    attr(ans, "details") <- ans.ret
    # print(ans)
    # sort by function value (DECREASING so best is last and follow.on gives natural ordering)
    if (ctrl$trace) { cat("Sort results\n") }
    if (ctrl$sort.result) { # sort by fvalues decreasing
	ord <- rev(order(as.numeric(ans$fvalues)))
	ans <- ans[ord, ]
    }
    if (ctrl$trace) { cat("DONE!\n") }

    return(ans)
} ## end of optimx
