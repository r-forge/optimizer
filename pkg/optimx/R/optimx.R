optimx <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=NULL,
            control=list(),
             ...) {
##### OPEN ISSUES: (any date order)
# 090729 -- ?? add minqa esp. BOBYQA
# 090729 -- ?? add Rcgminu and Rcgminb for large scale problems
# 090729 -- ?? simplify structure of answers -- avoid repetition and make easier to access
# 090531 -- ?? need SNewton back in to keep Ramsay and Co. on side
# 090531 -- ?? Use function code joint with funtest and funcheck in initial checks
# 090601 -- ?? Put optimx.dev into package to provide for local source packages etc.
# 090601 -- ?? Do we want hessevals to count Hessian evaluations?
# 090612 -- ?? There may be better choices for the tolerances in tests of equality for gradient / kkt tests.

##### IMPLEMENTED: (reverse date order)
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
#
#
#
# Output:
# ans = an object containing two sets of information:
# essential output = a data frame of converged parameters, function value at convergence, name of algorithm, 
#     convergence indicator, and function evaluation and gradient evaluation counts
# detailed output = this is a list with one component for each algorithm, and contains parameters, function values, convergence code, number of function and gradient evals, 
#	numerical gradient and hessian at convergence, eigenvalues of hessian, cpu time, other nformation returned by algorithms, and name of the algorithm.
# detailed output can be accessed via the attribute called `details'
#
#  Authors:  Ravi Varadhan & John Nash
#  Date:  February 17, 2008
#  Changes: Ravi Varadhan - Date: May 29, 2009
#
#################################################################
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
       all.methods=FALSE
    ) # for now turn off sorting until we can work out how to do it nicely
    
# Note that we do NOT want to check on the names, because we may introduce 
#    new names in the control lists of added methods
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
# However, we do want to substitute the appropriate information. 
# hessian control gets copied to kkt control
    if (!is.null(hessian)){
	control$kkt<-hessian
    }
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
       # print(onename)
       if (onename %in% nctrl) {
       #    cat("Existing\n") 
           ctrl[onename]<-control[onename]
       } else {
       #    cat("New\n")
           ctrl[onename]<-control[onename]
       }
    }
    if (!is.null(control$kkt)) { # default is to turn off kkt for large matrices
      if (is.null(gr)) {
         if (npar > 50) {
           ctrl$kkt=FALSE # too much work when large number of parameters
           if (ctrl$trace) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      } else {
         if (npar > 500) {
            ctrl$kkt=FALSE # too much work when large number of parameters
                    # even if we have a gradient
            if (ctrl$trace) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      }
    }

# Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) { have.bounds<-TRUE # set this for convenience
  } else { have.bounds <- FALSE }
  # print(have.bounds)

# Check parameters in bounds (090601: As yet not dealing with masks ??)
  #  bdmsk<-as.vector(bdmset[k, ])

  infeasible<-FALSE
  if (ctrl$trace > 0) cat("Function has ",npar," arguments\n")
 
  if (have.bounds) {
    # Expand bounds to vectors if needed
    if (length(lower)==1 ) lower <- rep(lower, npar)
    if (length(upper)==1 ) upper <- rep(upper, npar)
    bstate<-vector(mode="character", length=npar)
    print(par)
    print(lower)
    print(upper)
    for (i in 1:npar) {
       if ( (lower[i]<=par[i]) && (par[i]<=upper[i])) {
    	bstate[i]<-" In Bounds "
        } else { 
        #      if(bdmsk[i]!=0) {
        
                  infeasible<-TRUE
        #      }
              if (lower[i]>par[i]) {bstate[i]<-" Out of Bounds LOW" } else { bstate[i]<-" Out of Bounds HIGH " }
        } # end if in bounds
#        if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bdmsk[i],"   ",bstate,"\n")
         if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bstate,"\n")
    } # end of for loop over parameter vector elements
    if(infeasible) { 
#        print(bstate)
        stop("Infeasible point, no further tests") # for moment, just stop, but later need better exit mechanism
    } 
  } # end have.bounds

  # Check if function can be computed
  firsttry<-try(finit<-fn(par, ...) ) # Note: This incurs one EXTRA function evaluation because optimx is
                                      # a wrapper for other methods
  if (class(firsttry) == "try-error") {
    infeasible <- TRUE
    print(firsttry)
    stop("Cannot evaluate function at initial parameters")
    # Also check that it is returned as a scalar
    if (is.vector(finit) || is.list(finit) || is.matrix(finit) || is.array(finit) || ! is.numeric(finit) ) {
       stop("Function provided is not returning a scalar number")
    }
  }
  if (is.infinite(finit) || is.na(finit)) {
     stop("Function returned is infinite or NA (non-computable)")
  }

  # Check that we have the functions we need
  ipkgs <- installed.packages()
  if ("numDeriv" %in% ipkgs[,1]) require(numDeriv, quietly=TRUE) 
  else stop("Install package `numDeriv'", call.=FALSE)


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


  if(! is.null(hess)){ # check gradient
    hname <- deparse(substitute(hess))
    if (ctrl$trace) cat("Analytic hessian from function ",hname,"\n\n")
    hn <- hessian(func=fn, x=par,...) # ?? should we use dotdat
    ha <- hess(par, ...)
    # Now test for equality
    teps <- (.Machine$double.eps)^(1/3)
    if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
  } else if (ctrl$trace) cat("Analytic Hessian not made available.\n")

## >>> End of code common to funtest, funcheck and optimx

  # Set up the vector to return answers
  ans.ret <- vector(mode="list")
  # mode= is not strictly required. length defaults to 0. Sets up our answer vector.

  # List of methods in base or stats, namely those in optim(), nlm(), nlminb()
  bmeth <- c("BFGS", "CG", "Nelder-Mead", "SANN", "L-BFGS-B", "nlm", "nlminb")

  # List of methods in packages. 
  pmeth <- c("spg", "ucminf")

  allmeth <- c(bmeth, pmeth)
  # cat("start with allmeth=\n")
  # print(allmeth)

  nomaxmeth <- c("nlm", "nlminb", "ucminf") # methods not (yet) allowing maximize

  # Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) allmeth <- c("L-BFGS-B", "nlminb", "spg") 
  # cat("After check limits, allmeth=\n")
  # print(allmeth)


  if (ctrl$maximize) {
     for (meth in allmeth) {
         if (meth %in% nomaxmeth) { 
             if(ctrl$trace>0) cat("Method ",meth," cannot maximize -- dropping\n")
             allmeth<-allmeth[! (meth == allmeth)]
         }
     }
  }
  # cat("After check maximize, allmeth=\n")
  # print(allmeth)
 
  if (ctrl$all.methods) {
	method<-allmeth
	if (ctrl$trace>0) {
		cat("all.methods is TRUE -- Using all available methods\n")
		print(method)
	}
  }

  # Expand bounds to vectors if needed
  if (any(is.finite(lower)) & length(lower)==1) lower <- rep(lower, length(par))
  if (any(is.finite(upper)) & length(upper)==1) upper <- rep(upper, length(par))

  # Partial matching of method string allowed
  method <- unique(match.arg(method, allmeth, several.ok=TRUE) )
  nmeth <- length(method) # number of methods requested

  ## Check that methods are indeed available and loaded
  for (i in 1:nmeth) {
     cmeth <- method[i]
     if (ctrl$trace > 0) cat("Looking for method = ",cmeth,"\n")
     # look first in base/stats etc.
     if (! (cmeth %in% allmeth) ) {
         errmsg <- paste(cmeth," not found in any list of methods available")
         stop(errmsg, call.=FALSE)
     } # otherwise the method is available, and just needs to be loaded
  } # end check methods available
#  
  if (ctrl$trace>1) {
    cat("Methods:")
    print(method)
  }	
  if(any(method == "spg")) {
	if ("BB" %in% ipkgs[,1]) require(BB, quietly=TRUE)
	else  stop("Package `BB' Not installed", call.=FALSE)
	}
	
  if(any(method == "ucminf")) { 
	if ("ucminf" %in% ipkgs[,1]) require(ucminf, quietly=TRUE)
	else  stop("Package `ucminf' Not installed", call.=FALSE)
	}

  times <- rep(0, nmeth)
  ## figure out how many methods and reserve that many times to record.
  names(times) <- method
  ## And label these times appropriately with the method.
  j <- 0
  ## j will be the index to the answers, one per each method, starting parameter pair

  for (i in 1:nmeth) {
      ## loop over the methods
      meth <- method[i]
      ## extract the method name
      conv <- -1
      ## indicate that we have not yet converged
      if (ctrl$trace) cat("Method: ", meth, "\n") # display the method being used

      # Extract control information e.g., maxit
      # cat("Str of ctrl$maxit\n")
      # str(ctrl$maxit)
      if (! is.null(ctrl$maxit) ) {
        if (length(ctrl$maxit) == 1 ) {
          if(ctrl$trace>0) cat("setting meth.maxit to common value ",ctrl$maxit,"\n")
          meth.maxit <- ctrl$maxit # single value for all tries
        } else {
          if(ctrl$trace>0) cat("setting meth.maxit to choice ",i," = ",ctrl$maxit[i],"\n")
          if (length(ctrl$maxit) != nmeth) stop("ctrl$maxit has wrong number of elements", call.=FALSE)
          meth.maxit <- ctrl$maxit[i] # for the case with separate methods
        }
      } else { meth.maxit<-NULL }
      # create local control list
      mcontrol<-ctrl
      mcontrol[which(names(mcontrol) == "maxit")] <- meth.maxit
      mcontrol$follow.on<-NULL
      mcontrol$save.failures<-NULL
      mcontrol$sort.result<-NULL
      mcontrol$kkt<-NULL
      mcontrol$all.methods<-NULL
      if (meth != 'spg') mcontrol$maximize<-NULL # only spg has maximize
      if (ctrl$maximize) { # Try to sort out "maximize" control
         if (any(meth == nomaxmeth)) {
		stop("At the moment you must explicitly set up maximization for ", meth)
         }
         if (! is.null(ctrl$fnscale)) { 
 		stop("Mixing maximize and fnscale is dangerous. Please correct.")
         } 
      }


      if (meth=="Nelder-Mead" || meth == "BFGS" || meth == "L-BFGS-B" || meth == "CG" || meth == "SANN") {
        ## Take care of methods BFGS, L-BFGS-B, and CG  from optim()
        time <- system.time(ans <- try(optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, 
                method=meth, control=mcontrol, ...), silent=TRUE))[1]

        # The time is the index=1 element of the system.time for the process, 
        # which is a 'try()' of the regular optim() function
        # cat("Debug optim call: ans structure\n")
        # str(ans)
        if (class(ans) != "try-error") conv <- ans$convergence
	#      convergence: An integer code. '0' indicates successful convergence.
	#      ?? Put in other codes !!
        ans$fevals<-ans$counts[1]
        ans$gevals<-ans$counts[2]
        ans$counts<-NULL # and erase it now saved
        ans$niter<-NULL
      }   ## end if using BFGS etc.
## --------------------------------------------
      else if (meth == "nlminb") {
        ## Use portLib routine nlminb
        mcontrol$iter.max<-mcontrol$maxit # different name in this routine
        mcontrol$maxit<-NULL
        time <- system.time(ans <- try(nlminb(start=par, objective=fn, gradient=gr, lower=lower, 
		upper=upper, control=mcontrol,  ...), silent=TRUE))[1]
        ## Here we use nlminb rather than optim as our minimizer
        if (class(ans) != "try-error") conv <- ans$convergence
        # Translate output to common format
        ans$value<-ans$objective
        ans$objective<-NULL
	# if (conv == 0) names(ans)[2] <- "value" # We want value even if not conv. so this commented out
        ans$fevals<-ans$evaluations[1]
        ans$gevals<-ans$evaluations[2]
	ans$evaluations<-NULL # cleanup
        ans$niter<-ans$iterations
        ans$iterations<-NULL
      }  ## end if using nlminb
## --------------------------------------------
      else if (meth == "nlm") {
        ## Use stats package nlm
           if (!is.null(gr)) attr(fn, "gradient") <- gr(par, ...)
           if (!is.null(hess)) attr(fn, "hessian") <- hess(par, ...)
        cat("About to try ans <- try(nlm(start=par, objective=fn, ...)\n")
        time <- system.time(ans <- try(nlm(p=par, f=fn, ...), silent=TRUE))[1]
        print(ans)
        ## Here we use nlm rather than optim as our minimizer
        if (class(ans) != "try-error") {conv <- ans$code; if (conv == 1 | conv == 2 | conv == 3) conv <- 0}
        # Translate output to common format
        if (conv == 0) {names(ans)[1] <- "value"; names(ans)[2] <- "par"}
        ans$conv<-conv
        ans$fevals<-NA
        ans$gevals<-NA # ?? need to fix this somehow in nlm code
        ans$niter<-ans$iterations
        ans$iterations<-NULL
      }  ## end if using nlm
## --------------------------------------------
      else if (meth == "spg") {
        ## Use BB routine spg
        time <- system.time(ans <- try(spg(par=par, fn=fn, gr=gr, lower=lower, upper=upper,  
		control=list(maxit=1500, trace=FALSE), ...), silent=TRUE))[1]
        if (class(ans) != "try-error") conv <- ans$convergence
        ## Here we don't need to set names.
        ans$fevals<-ans$feval
        ans$feval<-NULL # to erase conflicting name
        ans$gevals<-NA # ??fixup needed
        ans$niter<-ans$iter
        ans$iter<-NULL
      }  ## end if using spg
## --------------------------------------------
      else if (meth == "ucminf") {
        ## Use ucminf routine
        mcontrol$meval<-mcontrol$maxit # different name in this routine
        mcontrol$maxit<-NULL
#        time <- system.time(ans <- try(ucminf(par=par, fn=fn, gr=gr,  control=mcontrol, ...), silent=TRUE))[1]
        time <- system.time(ans <- try(ucminf(par=par, fn=fn, gr=gr,  ...), silent=TRUE))[1]
        if (class(ans) != "try-error") {conv <- ans$convergence; if (conv == 1 | conv == 2 | conv == 4) conv <- 0}
# Termination criteria are tricky here!  How to determine successful convergence?
        ## Here we don't need to set names
        ans$fevals<-ans$info[4]
        ans$gevals<-ans$info[4] # calls fn and gr together
        ans$info<-NULL # to erase conflicting name
        ans$niter<-NULL
        if (ans$convergence == 1 || ans$convergence == 2 || ans$convergence == 4) {
		ans$conv <- 0
	} else {
		ans$conf <- ans$convergence
	}        
      }  ## end if using ucminf
## --------------------------------------------
## ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINDED METHOD: ", meth, sep='')
             stop(errmsg, call.=FALSE)
           }
## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      times[i] <- times[i] + time
      ## Accumulate time for a single method

      if ( ctrl$save.failures || (conv == 0) ){  ## We have a convergence, so save the solution
        ## Why do we not want to do anything with failures -- they may use even MORE time, and we may be "close"
        j <- j + 1 ## increment the counter for (successful) method/start case
        if (ctrl$trace) cat("Successful convergence! \n")  ## inform user we've succeeded

        # Testing final solution? Use numDeriv to get gradient and Hessian; compute Hessian eigenvalues
	ans$kkt1<-NA
	ans$kkt2<-NA
        # cat("Post processing -- ctrl$kkt = ",ctrl$kkt,"\n")
        if(ctrl$kkt) {
          ngatend <- if (is.null(gr)) grad(fn, ans$par) else gr(ans$par) # Gradient at solution
          nhatend <- if (is.null(gr)) hessian(fn, ans$par) else jacobian(gr,ans$par) # numerical hessian at "solution"
          # For bounds constraints, we need to "project" the gradient and Hessian
          bset<-sort(unique(c(which(ans$par<=lower), which(ans$par>=upper))))
          nbds<-length(bset) # number of elements nulled by bounds
  	  # Note: we assume that we are ON, not over boundary, but use <= and >=
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
        # Do we want more information??
        ans.ret[[j]] <- ans  ## save the answer. [[]] indexes the CONTENTS of the list
        ans.ret[[j]]$method <- method[i] ## and we tag on the method with the $ linker
      }  ## end post-processing of successful solution

	if (ctrl$follow.on) par <- ans$par
    } ## end loop over method (index i)
    attr(ans.ret, "CPU times (s)") <- times ## save the accumulated times 
    pars <- lapply(ans.ret, function(x) x$par)
    vals <- lapply(ans.ret, function(x) x$value)
    meths <- lapply(ans.ret, function(x) x$method)
    fevals<- lapply(ans.ret, function(x) x$fevals)
    gevals<- lapply(ans.ret, function(x) x$gevals)
    nitns <-  lapply(ans.ret, function(x) x$niter)
    convcode<- lapply(ans.ret, function(x) x$conv)
    kkt1<- lapply(ans.ret, function(x) x$kkt1)
    kkt2<- lapply(ans.ret, function(x) x$kkt2)

    ans <- data.frame(cbind(par=pars, fvalues=vals, method=meths, fns=fevals, grs=gevals, 
	itns=nitns, conv=convcode, KKT1=kkt1, KKT2=kkt2))
    attr(ans, "details") <- ans.ret

    # sort by function value (DECREASING so best is last and follow.on gives natural ordering)
    if (ctrl$sort.result) { # sort by fvalues decreasing
	ord <- rev(order(as.numeric(ans$fvalues)))
	ans <- ans[ord, ]
    }
    return(ans)
} ## end of optimx


