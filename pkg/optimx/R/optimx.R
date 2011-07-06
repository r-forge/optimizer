optimx <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, bdmsk=NULL, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {
##### OPEN ISSUES: (any date order)
# 110703 -- ?? check maximize works in examples
# 110701 -- ?? move scalecheck into bmchk package or into a scale check package??
# 110628 -- ?? check how tfn is built for nlm() with ugr and uhess when gr or hess are NULL
# 110628 -- ?? put gHgenb into starttests and replace hessian()
# 110628 -- ?? tests needed of hess with maximize
# 110628 -- ?? tests needed with bounds and masks active
# 110628 -- ?? return bdmsk[]. Where? In "details" or in main answer?
# 110625 -- ?? dealing with bounds and masks with gHgen and kktc -- gHgenb
# 110625 -- ?? dealing with bounds and masks in numDeriv and simple numerical derivs
# 110624 -- ?? align starttests with funcheck package to avoid maintaining double code
# 110615 -- ?? should we put more in ufn e.g., function eval count etc.
# 110615 -- ?? should we put in ugr like we have ufn
# 091220 -- ?? kkt1 tolerance. Do we want to let this be an input control?
# 090612 -- ?? (Similar to above) Better choices for the tolerances in tests of equality for gradient / kkt tests?
# 091018 -- ?? masks? generally in optimx
# 090929 -- ?? control structure could use simplification 
# 090729 -- ?? simplify structure of answers -- avoid repetition and make easier to access
# 090531 -- ?? like SNewton in to assist other workers
# 090601 -- ?? Do we want hessevals to count Hessian evaluations?

##### IMPLEMENTED: (reverse date order)
# 110701 -- combined fnuser has fnuser$fn, fnuser$gr, fnuser$hess, of which 
#           3rd or 2nd and 3rd may be NULL. Still do NOT have right sense of how
#           to put in appropriate numerical approximations to these. ??
# 110627 -- hess() check does not look at whether we maximize -- fixed
# 110624 -- put in gHgen and kktc from kktc package
# 110615 -- sorting answers when maximize=TRUE
# 110614 -- put in myfn.R as ufn to watch for inadmissible / non-computable function
# 110212 -- & to && and | to || for controls
# 110212 -- Hessian changed from NULL to FALSE default (??explain)
# 100328 -- check if maximize works for Powell (minqa) routines -- 100415
# 100329 -- make sure spg fixed on CRAN -- 100415
# 100329 -- maximize tested for all but minqa, though spg needs fixup on CRAN. ??
# 100215 -- Add newuoa and uobyqa to test minqa
# 100212 -- setting scaletol in controls
# 091018 -- new control "starttest" so we can skip them DONE 091220
# 091218 -- scaling checks (before computing rhobeg) and warning
# 091215 -- mcontrol omission in nlm() call fixed
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
#                 The choice of BFGS second is open to discussion. JN
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
#      sort.result=TRUE, that is, we sort the results in decreasing order of the final function value
#      kkt=TRUE to run Kuhn, Karush, Tucker tests of results unless problem large
#      kkttol=0.001 (was .Machine$double.eps^(1/4)) Default value to check for small gradient and negative
#               Hessian eigenvalues
#      kkt2tol=1E-6 (default WAS 10* default for kkttol) Tolerance for eigenvalue ratio in KKT test of 
#               positive definite Hessian
#      all.methods=FALSE By default we do NOT run all methods
#      starttests=TRUE  By default we run tests of the function and parameters: feasibility relative to
#               bounds, analytic vs numerical gradient, scaling tests) before we try optimization methods
#      dowarn=TRUE By default we leave warnings generated by optimx.
#      badval=(0.5)*.Machine$double.xmax The value to set for the function value when try(fn()) fails.
#      scaletol=3 To check parameters or their bounds we take logs of absolute values and find the range 
#               i.e., max - min. This should not exceed scaletol. A value of 3 gives magnitudes between 
#               1 and 20 approximately.  
#
# Output:
# ans = an object containing two sets of information:
# essential output = a data frame of converged parameters, function value at convergence, 
#    name of algorithm, convergence indicator, and function evaluation and 
#    gradient evaluation counts
# detailed output = this is a list with one component for each algorithm, and contains parameters, 
#    function values, convergence code, number of function and gradient evals, numerical gradient 
#    and hessian at convergence, eigenvalues of that hessian approximation, cpu time, and other 
#    information returned by algorithms, and name of the algorithm.
# detailed output can be accessed via the attribute called `details'
#
#  Authors:  Ravi Varadhan & John Nash
#  Date:  February 17, 2008
#  Changes: Ravi Varadhan - Date: May 29, 2009, John Nash - Latest: July 2, 2011
#
#################################################################
nullgr<-is.null(gr) # save these as we redefine functions so NOT null later
nullhess<-is.null(hess)
pbscale<-scalecheck(par, lower = lower, upper = upper, bdmsk=NULL, dowarn = TRUE)

# Get real name of function to be minimized
   fname<-deparse(substitute(fn))
   if (is.null(control$trace)) control$trace<-0 # to ensure trace set
   if (control$trace>0) {
      cat("fn is ",fname,"\n")
   }
## Code more or less common to funtest, funcheck and optimx <<<
# Check parameters are in right form
   if(!is.null(dim(par))) stop("Parameter should be a vector, not a matrix!")
   if (! is.vector(par) ) {
      stop("The parameters are NOT in a vector")
   }
   npar<-length(par)
   if (is.null(bdmsk) ) bdmsk<-rep(1,npar) # set masks for free parameters
# Set control defaults. Do we want to save.failures?
   ctrl <- list(
   follow.on=FALSE, 
   save.failures=TRUE,
   eps=1.e-07,
   trace=0,
   sort.result=TRUE,
   kkt=TRUE,
   all.methods=FALSE,
   starttests=TRUE,
   maximize=FALSE,
   dowarn=TRUE, 
   kkttol=0.001,
   kkt2tol=1.0E-6,
   badval=(0.5)*.Machine$double.xmax,
   scaletol=3,
   usenumDeriv=FALSE
   ) # for now turn off sorting until we can work out how to do it nicely
# Note that we do NOT want to check on the names, because we may introduce 
#    new names in the control lists of added methods. That is, we do NOT do
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
# However, we do want to substitute the appropriate information. 
# removed copy of hessian to control$kkt
   ncontrol <- names(control)
   nctrl <- names(ctrl)
   for (onename in ncontrol) {
      if (!(onename %in% nctrl)) {
         if (ctrl$trace>2) cat("control ",onename," is not in default set\n")
      }
      ctrl[onename]<-control[onename]
   }
   if (is.null(control$kkt)) { # turn off kkt for large matrices
      ctrl$kkt<-TRUE # default it to compute KKT tests
      if (is.null(gr)) { # no analytic gradient
         if (npar > 50) {
            ctrl$kkt=FALSE # too much work when large number of parameters
            if (ctrl$trace>0) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      } else {
         if (npar > 500) {
            ctrl$kkt=FALSE # too much work when large number of parameters, even with analytic gradient
            if (ctrl$trace>0) cat("gr NULL, npar > 50, kkt set FALSE\n")
         }
      }
   } else { # kkt is set
      if (control$kkt) {
         if (nullgr) {
            if (npar > 50) {
               if ((ctrl$trace>0) && ctrl$dowarn) warning("Computing hessian for gr NULL, npar > 50, can be slow\n")
            }
         } else {
            if (npar > 500) {
               if ((ctrl$trace>0) && ctrl$dowarn) warning("Computing hessian with gr code, npar > 500, can be slow\n")
            }
         }
      }
   }
   # copy some controls to top level
   trace<-ctrl$trace
   maximize<-ctrl$maximize
   dowarn<-ctrl$dowarn
   eps<-ctrl$eps
   # end move ctrl elements to top level
   if ((trace>0) && maximize ) cat("MAXIMIZING\n")
   #############################################
    fnuser<-list(fn=fn, gr=gr, hess=hess) # define the user function
   #############################################
# Do we have bounds? ?? move to bmchk!!
   if (any(is.finite(c(lower, upper)))) { have.bounds<-TRUE # set this for convenience
   } else { have.bounds <- FALSE }
   if (ctrl$starttests) {
      # Check parameters in bounds (090601: As yet not dealing with masks ??)
      infeasible<-FALSE
      if (ctrl$trace > 0) cat("Function has ",npar," arguments\n")
      if (have.bounds) {
         # Expand bounds to vectors if needed
         if (length(lower)==1) lower <- rep(lower, npar)
         if (length(upper)==1) upper <- rep(upper, npar)
         if (length(lower)!=npar) stop("lower bounds vector of wrong length")
         if (length(upper)!=npar) stop("upper bounds vector of wrong length")
         bstate<-vector(mode="character", length=npar)
         for (i in 1:npar) {
            if ( (lower[i]<=par[i]) && (par[i]<=upper[i]) ) {
               bstate[i]<-" In Bounds "
            } else { 
               if (bdmsk[i]!=0) {
                  infeasible<-TRUE
               }
               if (lower[i]>par[i]) {bstate[i]<-" Out of Bounds LOW " } else { bstate[i]<-" Out of Bounds HIGH " }
            } # end if in bounds
            if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bstate[i],"\n")
         } # end of for loop over parameter vector elements
         if (infeasible) { 
            stop("Infeasible point, no further tests")
         } 
      } # end have.bounds
      # Check if function can be computed  
      finit<-ufn(par, fnuser, maximize, ...) # 110626
      # Note: This incurs one EXTRA function evaluation because optimx wraps other methods
      if (attr(finit, "inadmissible")) {
         infeasible <- TRUE
         stop("Cannot evaluate function at initial parameters")
	##?? eventually want to exit WITHOUT stopping to maintain program control
        ## but WITH a proper exit code
      }
      # Also check that it is returned as a scalar
      if (is.vector(finit)) {
         if ( length(finit)==1 )  stop("!(is.vector(finit) && (length(finit)==1) )")
      }
      if (is.list(finit)) stop("Initial function is a list!")
      if (is.matrix(finit)) stop("Initial function is a matrix!")
      if ( is.array(finit) ) stop("Initial function is an array!")
      if (  ! (is.numeric(finit)) ) stop("Initial function is not numeric")
      if (is.infinite(finit) || is.na(finit)) {
         stop("Function returned is infinite or NA (non-computable)")
      }
   } # end have starttests -- computable function
   if (ctrl$trace>0) cat("Function at initial point=",finit,"\n")
   # Check that we have the functions we need
   ipkgs <- installed.packages()
   if ("numDeriv" %in% ipkgs[,1]) require(numDeriv, quietly=FALSE) 
   else stop("Install package `numDeriv'", call.=FALSE)
   # gradient etc. start tests -- move to package funtest!!!
   if(ctrl$starttests) {
      if(! nullgr){ # check gradient
         gname <- deparse(substitute(gr)) # ??? but embedded in ugr
         if (ctrl$trace>1) cat("Analytic gradient uses function ",gname,"\n")
         fval<-finit # already evaluated
#         fval <- ufn(par, fn, ctrl$maximize, ...)
         if (ctrl$trace>1) cat("Compute analytic gradient\n")
         ga <- ugr(par, fnuser, ctrl$maximize, ... ) # analytic gradient
         if (ctrl$trace>1) print(ga)
         if (ctrl$trace>1) cat("Compute numeric gradient\n")
         gn <- grad(func=ufn, x=par, fnuser=fnuser, maximize=ctrl$maximize, ...) # numerically approximated gradient
         if (ctrl$trace>1) print(gn)
         # Now test for equality (090612: ?? There may be better choices for the tolerances.
         teps <- (.Machine$double.eps)^(1/3)
         if (max(abs(gn-ga))/(1 + abs(fval)) >= teps) 
               stop("Gradient function might be wrong - check it! \n", call.=FALSE)
      } else if (ctrl$trace>0) cat("Analytic gradient not made available.\n")
      if(! nullhess){ # check Hessian (?? what about using gHgenb or funtest??)
         hname <- deparse(substitute(hess))
         if (ctrl$trace>0) cat("Analytic hessian from function ",hname,"\n\n")
         ha <- uhess(par, fnuser, ctrl$maximize, ... ) # analytic hessian 
         if (attr(ha,"inadmissible")) {
		stop("Cannot evaluate starting analytic Hessian")
		# want to be able to keep running??
         }
         if (nullgr) {
            hn <- hessian(func=ufn, x=par, maximize=ctrl$maximize, fnuser=fnuser, ...) 
         } else {
            hn <- jacobian(func=ugr, x=par, maximize=ctrl$maximize, fnuser=fnuser, ...)
            hn <- 0.5*(hn+t(hn)) # force symmetry
         }
         # Now test for equality ?? again have to consider tolerance
         teps <- (.Machine$double.eps)^(1/3)
         if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
      } else if (ctrl$trace>0) cat("Analytic Hessian not made available.\n")
   } # end gradient/hessian tests
## >>> End of code common to funtest, funcheck and optimx, with mods for local needs
# Set up the vector to return answers
   ans.ret <- vector(mode="list")
   # mode= is not strictly required. length defaults to 0. This sets up our answer vector.
# List of methods in base or stats, namely those in optim(), nlm(), nlminb()
   bmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb")
# SANN has no termination for optimality, only a maxit count for
#    the maximum number of function evaluations; remove DEoptim for now -- not useful 
#    for smooth functions. Code left in as example for those who may need it.
# List of methods in packages. 
# uobyqa removed 110114 because of some crashes that did not seem resolvable. newuoa should be more efficient anyway.  
   pmeth <- c("spg", "ucminf", "Rcgmin", "Rvmmin", "bobyqa", "newuoa")
   allmeth <- c(bmeth, pmeth)
   # Restrict list of methods if we have bounds
   if (any(is.finite(c(lower, upper)))) allmeth <- c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa") 
   if (ctrl$all.methods) { # Changes method vector!
      method<-allmeth
      if (ctrl$trace>0) {
         cat("all.methods is TRUE -- Using all available methods\n")
         print(method)
      }
   } 
   # Methods that can allow maximization
   if (ctrl$maximize) {
      if (ctrl$trace>3) cat("remove nlm as we are maximizing")
      method<-method[which(method != "nlm")]
   }
   # Partial matching of method string allowed
   # avoid duplicates here
   # 2011-1-17 JN: to set L-BFGS-B
   method <- try(unique(match.arg(method, allmeth, several.ok=TRUE) ),silent=TRUE)
   if (class(method)=="try-error") {
      warning("optimx: No match to available methods")
      method<-NULL
      nmeth<-0
   } else {
      nmeth <- length(method) # number of methods requested
   } 
   # JN 2011-1-17 fix for default when there are bounds
   if ((nmeth==0) && have.bounds) {
      method="L-BFGS-B"
      warning("Default method when bounds specified is L-BFGS-B to match optim()")
      nmeth<-1
   }
   if (nmeth==0) stop("No suitable method")
   if (ctrl$trace>1) {
      cat("Methods to be used:")
      print(method)
   }

   ## Check that methods are indeed available and loaded
   for (i in 1:nmeth) {
      cmeth <- method[i]
      if (ctrl$trace > 1) cat("Looking for method = ",cmeth,"\n")
      if (! (cmeth %in% allmeth) ) {
         errmsg <- paste(cmeth," not found in any list of methods available")
         stop(errmsg, call.=FALSE)
      } # otherwise the method is available, and just needs to be loaded
   } # end check methods available
   # restrict methods if we have masks
   if (any( bdmsk == 0 ) ){
      mmeth<-c("Rcgmin", "Rvmmin")
      if (ctrl$trace > 0) cat("Masks in effect. Restrict methods")
      method <- try(unique(match.arg(method, mmeth, several.ok=TRUE) ),silent=TRUE)
      if (class(method)=="try-error") {
         warning("optimx: No match between requested and mask-capable methods")
         method<-mmeth # THIS MAY BE UNWANTED!
         nmeth <- length(method) 
      } else {
         nmeth <- length(method) # number of methods requested
      }
   }
# Now make sure methods loaded
   if(any(method == "spg")) {
      if ("BB" %in% ipkgs[,1]) require(BB, quietly=FALSE)
      else  stop("Package `BB' Not installed", call.=FALSE)
   }
   if(any(method == "ucminf")) { 
      if ("ucminf" %in% ipkgs[,1]) require(ucminf, quietly=FALSE)
      else  stop("Package `ucminf' Not installed", call.=FALSE)
   }
   if(any(method == "Rcgmin")) { 
      if ("Rcgmin" %in% ipkgs[,1]) require(Rcgmin, quietly=FALSE)
      else  stop("Package `Rcgmin' Not installed", call.=FALSE)
   }
   if(any(method == "Rvmmin")) { 
      if ("Rvmmin" %in% ipkgs[,1]) require(Rvmmin, quietly=FALSE)
      else  stop("Package `Rvmmin' Not installed", call.=FALSE)
   }
   if(any(method == "bobyqa")) { 
      if ("minqa" %in% ipkgs[,1]) require(minqa, quietly=FALSE)
      else  stop("Package `minqa' (for bobyqa) Not installed", call.=FALSE)
   }
   if(any(method == "bobyqa")) { 
      if ("minqa" %in% ipkgs[,1]) require(minqa, quietly=FALSE)
      else  stop("Package `minqa' (for bobyqa) Not installed", call.=FALSE)
   }
   if(any(method == "uobyqa")) { 
      if ("minqa" %in% ipkgs[,1]) require(minqa, quietly=FALSE)
      else  stop("Package `minqa' (for uobyqa) Not installed", call.=FALSE)
   }
   if(any(method == "newuoa")) { 
      if ("minqa" %in% ipkgs[,1]) require(minqa, quietly=FALSE)
      else  stop("Package `minqa' (for newuoa) Not installed", call.=FALSE)
   }
#  if(any(method == "DEoptim")) { # Code removed as DEoptim not part of current set of methods
#     if ("DEoptim" %in% ipkgs[,1]) require(DEoptim, quietly=FALSE)
#     else  stop("Package `DEoptim' Not installed", call.=FALSE)
#  }
# Scaling check  091219
   if (ctrl$starttests) {
   srat<-scalecheck(par, lower, upper,ctrl$dowarn)
   sratv<-c(srat$lpratio, srat$lbratio)
   if(max(sratv,na.rm=TRUE) > ctrl$scaletol) { 
      warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
      if (ctrl$dowarn) warning(warnstr)
   }
        if(ctrl$trace>0) {
      cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
   }
    }
# end scaling check
# Run methods
   times <- rep(0, nmeth)  # figure out how many methods and reserve that many times to record.
   names(times) <- method  # And label these times appropriately with the method.
   j <- 0  # j will be the index to the answers, one per each method, starting parameter pair
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
         if (ctrl$follow.on) cat("Do ",ctrl$maxit," steps of ",meth,"\n")
      }
      if (ctrl$trace>0) cat("Method: ", meth, "\n") # display the method being used
      # Extract control information e.g., trace
      # 20100215: Note that maxit needs to be defined other than 0 e.g., for ucminf
      # create local control list for a single method -- this is one of the key issues for optimx
      mcontrol<-ctrl
      mcontrol$follow.on<-NULL # And make sure that controls NOT in method are deleted (nulled)
      mcontrol$save.failures<-NULL
      mcontrol$sort.result<-NULL
      mcontrol$kkt<-NULL
      mcontrol$starttests<-NULL
      mcontrol$all.methods<-NULL
      mcontrol$dowarn<-NULL
      mcontrol$kkttol<-NULL
      mcontrol$kkt2tol<-NULL
      mcontrol$maximize<-NULL # Even if method DOES have it
      mcontrol$badval<-NULL
      mcontrol$eps<-NULL
      mcontrol$scaletol<-NULL # not used in any methods -- it is here for the scale check of parameters and bounds above
# Methods from optim()
# ---  UNDEFINED METHOD ---
      if (! (meth %in% allmeth )) {
         errmsg<-paste("UNDEFINED METHOD: ", meth, sep='')
         stop(errmsg, call.=FALSE)
      }
## --------------------------------------------
      if (meth=="Nelder-Mead" || meth == "BFGS" || meth == "L-BFGS-B" || meth == "CG") { # no SANN
         if (ctrl$trace>2) cat("Found an optim method\n")
         # Take care of methods   from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
         mcontrol$usenumDeriv<-NULL # not used in optim()
         if (nullgr) {
             # cat("NULL gr\n")
             myfn<-list(fn=fn, gr=NULL, hess=NULL)
             time <- system.time(ans <- try(optim(par=par, fn=ufn, gr=NULL, 
                     lower=lower, upper=upper, method=meth, 
                     control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
             myfn<-list(fn=fn, gr=gr, hess=NULL)
              time <- system.time(ans <- try(optim(par=par, fn=ufn, gr=ugr, 
                     lower=lower, upper=upper, method=meth, 
                     control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         }
         # The time is the index=1 element of the system.time for the process, 
         # which is a 'try()' of the regular optim() function
         if (class(ans)[1] != "try-error") { 
            ans$conv <- ans$convergence
         #  convergence: An integer code. '0' indicates successful convergence.
         #  ans$value alread defined OK in this case
            ans$fevals<-ans$counts[1] # save function and gradient count information
            ans$gevals<-ans$counts[2]
            ans$counts<-NULL # and erase the counts element now data is saved
            ans$niter<-NULL # not used, so delete
         } else { # bad result -- What to do?
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$conv<-9999 # failed in run
            if (ctrl$trace>0) cat("optim function evaluation failure\n")
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
         }
         if (ctrl$maximize) { # want to return the maximum with correct sign
            cat("optim:\n")
            print(ans)
            ans$value= -ans$value
         }
      }   # end if using optim() methods
## --------------------------------------------
      else 
      if (meth == "nlminb") {
         # Here we use portLib routine nlminb rather than optim as our minimizer
         mcontrol$iter.max<-mcontrol$maxit # different name for iteration limit in this routine
         mcontrol$maxit<-NULL
         mcontrol$usenumDeriv<-NULL
         mcontrol$abs.tol<-0 # To fix issues when minimum is less than 0. 20100711
         if( is.null(mcontrol$trace) || is.na(mcontrol$trace) || mcontrol$trace == 0) { 
            mcontrol$trace = 0
         } else { 
            mcontrol$trace = 1 # this is EVERY iteration. nlminb trace is freq of reporting.
         }

         if (nullgr) {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(nlminb(start=par, objective=ufn, gradient=NULL, lower=lower, 
              upper=upper, control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(nlminb(start=par, objective=ufn, gradient=ugr,  lower=lower,
              upper=upper, control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         }
         if (class(ans)[1] != "try-error") {
            ans$conv <- ans$convergence
            # Translate output to common format and names
            ans$value<-ans$objective
            ans$objective<-NULL
            ans$fevals<-ans$evaluations[1]
            ans$gevals<-ans$evaluations[2]
            ans$evaluations<-NULL # cleanup
            ans$niter<-ans$iterations
            ans$iterations<-NULL
         } else { # bad result -- What to do?
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$conv<-9999 # failed in run
            if (ctrl$trace>0) cat("nlminb function evaluation failure\n")
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
            if (ctrl$trace) {
               cat("maximize using nlminb:\n")
               print(ans)
            }
         }
      }  ## end if using nlminb
## --------------------------------------------
#      else 
      if (meth == "nlm") { # Use stats package nlm routine
         tufn <- fn
         attr(tufn, "gradient") <- gr
         attr(tufn, "hessian") <- hess
         if (dowarn) warning("Method nlm function evaluation NOT safeguarded")
#         tufn <- ufn
#         if (nullgr) attr(tufn, "gradient")<-NULL else attr(tufn, "gradient") <- ugr
#         if (nullhess) attr(tufn, "hessian")<-NULL else attr(tufn, "hessian") <- uhess
 
         ## 091215 added control for iteration limit
         if (! is.null(mcontrol$maxit)) { iterlim<-mcontrol$maxit } 
         else { 
            iterlim = 100 
            mcontrol$maxit<-NULL # and remove it for this method
         }
         if (! is.null(mcontrol$trace)) { plevel<-0 }
         else {
            plevel<-min(2,mcontrol$trace) # 110702 max is 2
            mcontrol$trace<-NULL
         }
         nampar<-names(par) # save names 110508
         time <- system.time(ans <- try(nlm(f=tufn, p=par, iterlim=iterlim, 
                     print.level=plevel,...), silent=TRUE))[1]
         if (class(ans)[1] != "try-error") {
            ans$conv <- ans$code
            if (ans$conv == 1 || ans$conv == 2 || ans$conv == 3) ans$conv <- 0
            if (ans$conv == 4) ans$conv <- 1
            # Translate output to common format
            ans$value<-ans$minimum
            ans$par<-ans$estimate
            ans$estimate<-NULL
            ans$minimum<-NULL
            ans$fevals<-NA
            ans$gevals<-NA # ?? need to fix this somehow in nlm code
            ans$niter<-ans$iterations
            ans$iterations<-NULL
            names(ans$par)<-nampar # restore names 110508
         } else {
            if (ctrl$trace > 0) cat("nlm failed for this problem\n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      } # end if using nlm
## --------------------------------------------
#      else 
      if (meth == "spg") { # Use BB package routine spg as minimizer
         mcontrol$usenumDeriv<-NULL
         mcontrol$maximize<-NULL # Use external maximization approach
#        No! Use original functions and maximize in controls.
#         time <- system.time(ans <- try(spg(par=par, fn=fn, gr=gr, lower=lower, upper=upper,  
#              control=mcontrol, ...), silent=TRUE))[1]
         if (nullgr) {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(spg(par=par, fn=ufn, gr=NULL, lower=lower, upper=upper,  
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
            myfn<-list(fn=fn, gr=gr, hess=NULL)
            time <- system.time(ans <- try(spg(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper,  
              fnuser=myfn, maximize=ctrl$maximize, ..., control=mcontrol), silent=TRUE))[1]
         }
         if (class(ans)[1] != "try-error") { 
         ans$conv <- ans$convergence
         ans$fevals<-ans$feval
         ans$feval<-NULL # to erase conflicting name
         ans$gevals<-NA # ??fixup needed
         ans$niter<-ans$iter
         ans$iter<-NULL
      } else { # spg failed
         if (ctrl$trace > 0) cat("spg failed for this problem\n")
         ans<-list(fevals=NA) # ans not yet defined, so set as list
         ans$value= ctrl$badval
         ans$par<-rep(NA,npar)
         ans$conv<-9999 # failed in run
         ans$gevals<-NA 
         ans$niter<-NA
      }
      if (ctrl$maximize) {
         ans$value= -ans$value
      }
   }  # end if using spg
## --------------------------------------------
#      else 
      if (meth == "ucminf") {
         ## Use ucminf routine
         if (is.null(ctrl$maxit)) mcontrol$maxit<-500 # ensure there is a default value
         # Change 20100415 to avoid setting ctrl values when all.methods
         mcontrol$maxeval<-mcontrol$maxit # Note it is just function evals for ucminf
         mcontrol$maxit<-NULL
         mcontrol$usenumDeriv<-NULL
         if (nullgr) {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(ucminf(par=par, fn=ufn, gr=NULL, 
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
            myfn<-list(fn=fn, gr=gr, hess=NULL)
            time <- system.time(ans <- try(ucminf(par=par, fn=ufn, gr=ugr, 
              fnuser=myfn, maximize=ctrl$maximize, ..., control=mcontrol), silent=TRUE))[1]
         }
##         time <- system.time(ans <- try(ucminf(par=par, fn=fn, gr=gr,  control=mcontrol, ...), silent=TRUE))[1]
         if (class(ans)[1] != "try-error") {
            ans$conv <- ans$convergence
            # From ucminf documentation:  
            # convergence = 1 Stopped by small gradient (grtol).
            #               2 Stopped by small step (xtol).
            #               3 Stopped by function evaluation limit (maxeval).
            #               4 Stopped by zero step from line search
            #              -2 Computation did not start: length(par) = 0.
            #              -4 Computation did not start: stepmax is too small.
            #              -5 Computation did not start: grtol or xtol <= 0.
            #              -6 Computation did not start: maxeval <= 0.
            #              -7 Computation did not start: given Hessian not pos. definite.
            #      message: String with reason of termination.
         if (ans$conv == 1 || ans$conv == 2 || ans$conv == 4) {
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
         } else { # ucminf failed
            if (ctrl$trace > 0) cat("ucminf failed for this problem\n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using ucminf
## --------------------------------------------
###### progress point #########
#      else 
      if (meth == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
         bdmsk<-rep(1,npar)
         mcontrol$trace<-NULL
         mcontrol$maximize<-NULL
         if (ctrl$trace>0) mcontrol$trace<-1
# use original functions -- protected inside (NOT protected inside)
#         time <- system.time(ans <- try(Rcgmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, bdmsk=bdmsk,
#              control=mcontrol, ...), silent=TRUE))[1]
         if (nullgr) {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(Rcgmin(par=par, fn=ufn, gr=NULL, lower=lower, upper=upper, bdmsk=bdmsk,
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
            myfn<-list(fn=fn, gr=gr, hess=NULL)
            time <- system.time(ans <- try(Rcgmin(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper, bdmsk=bdmsk,
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         }
         if (class(ans)[1] != "try-error") {
            ans$conv <- ans$convergence
            ans$fevals<-ans$counts[1]
            ans$gevals<-ans$counts[2]
            # ans$value<-ans$value 
         } else {
            if (ctrl$trace>0) cat("Rcgmin failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using Rcgmin
## --------------------------------------------
###### progress point #########
#      else 
      if (meth == "Rvmmin") { # Use Rvmmin routine (ignoring masks)
         bdmsk<-rep(1,npar)
         mcontrol$trace<-NULL
         mcontrol$maximize<-NULL
         if (ctrl$trace>0) mcontrol$trace<-1
         mcontrol$maximize<-NULL # negation built into ufn
#        Use original functions (ufn ugr uhess inside)
         if (nullgr) {
            myfn<-list(fn=fn, gr=NULL, hess=NULL)
            time <- system.time(ans <- try(Rvmmin(par=par, fn=ufn, gr=NULL, lower=lower, upper=upper, bdmsk=bdmsk,
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         } else {
            myfn<-list(fn=fn, gr=gr, hess=NULL)
            time <- system.time(ans <- try(Rvmmin(par=par, fn=ufn, gr=ugr, lower=lower, upper=upper, bdmsk=bdmsk,
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         }
         if ((class(ans)[1] != "try-error") && (ans$convergence==0)) {
            ans$conv <- ans$convergence
            ans$fevals<-ans$counts[1]
            ans$gevals<-ans$counts[2]
         } else {
            if (ctrl$trace>0) cat("Rvmmin failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using Rvmmin
## --------------------------------------------
#      else 
      if (meth == "bobyqa") {# Use bobyqa routine from minqa package
         if (! is.null(mcontrol$maxit)) { 
            mcontrol$maxfun<-mcontrol$maxit
         } else {
            mcontrol$maxfun<-5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
         }
         mcontrol$iprint<-min(mcontrol$trace, 3)
         mcontrol$maxit<-NULL # and null out control that is NOT used
         mcontrol$trace<-NULL
         mcontrol$usenumDeriv<-NULL
         nampar<-names(par) # save names 110508
         myfn<-list(fn=fn, gr=NULL, hess=NULL)
         time <- system.time(ans <- try(bobyqa(par=par, fn=ufn, lower=lower, upper=upper, 
               control=mcontrol, fnuser=myfn, maximize=ctrl$maximize, ...), silent=TRUE))[1]
         if (class(ans)[1] != "try-error") {
            ans$conv <- 0
            if (ans$feval > mcontrol$maxfun) {
               ans$conv <- 1 # too many evaluations
            }
            ans$fevals<-ans$feval
            ans$gevals<-NA
            ans$value<-ans$fval 
            names(ans$par)<-nampar # restore names 110508
         } else {
            if (ctrl$trace>0) cat("bobyqa failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using bobyqa
## --------------------------------------------
#      else 
      if (meth == "uobyqa") { # Use uobyqa routine from minqa package
         if (! is.null(mcontrol$maxit)) { 
            mcontrol$maxfun<-mcontrol$maxit
         } else {
            mcontrol$maxfun<-5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
         }
         mcontrol$maxit<-NULL # and null out control that is NOT used
         mcontrol$iprint<-min(mcontrol$trace, 3)
         mcontrol$trace<-NULL
         mcontrol$usenumDeriv<-NULL
         nampar<-names(par) # save names 110508
         myfn<-list(fn=fn, gr=NULL, hess=NULL)
#         time <- system.time(ans <- try(uobyqa(par=par, fn=ufn, control=mcontrol,...), silent=TRUE))[1]
         time <- system.time(ans <- try(uobyqa(par=par, fn=ufn, 
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize,...), silent=TRUE))[1]
         if (class(ans)[1] != "try-error") {
            ans$conv <- 0
            if (ans$feval > mcontrol$maxfun) {
               ans$conv <- 1 # too many evaluations
            }
            ans$fevals<-ans$feval
            ans$gevals<-NA
            ans$value<-ans$fval 
            names(ans$par)<-nampar # restore names 110508
         } else {
            if (ctrl$trace>0) cat("uobyqa failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using uobyqa
## --------------------------------------------
#      else 
      if (meth == "newuoa") {# Use newuoa routine from minqa package
         if (! is.null(mcontrol$maxit)) { 
            mcontrol$maxfun<-mcontrol$maxit
         } else {
            mcontrol$maxfun<-5000*round(sqrt(npar+1)) # ?? default at 100215, but should it be changed?!!
         }
         mcontrol$maxit<-NULL # and null out control that is NOT used
         mcontrol$iprint<-min(mcontrol$trace, 3)
         mcontrol$trace<-NULL
         mcontrol$usenumDeriv<-NULL
         nampar<-names(par) # save names 110508
         myfn<-list(fn=fn, gr=NULL, hess=NULL)
         time <- system.time(ans <- try(newuoa(par=par, fn=ufn, 
              control=mcontrol, fnuser=myfn, maximize=ctrl$maximize,...), silent=TRUE))[1]
         if (class(ans)[1] != "try-error") {
            ans$conv <- 0
            if (ans$feval > mcontrol$maxfun) {
               ans$conv <- 1 # too many evaluations
            }
            ans$fevals<-ans$feval
            ans$gevals<-NA
            ans$value<-ans$fval 
            names(ans$par)<-nampar # restore names 110508
         } else {
            if (ctrl$trace>0) cat("newuoa failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value= ctrl$badval
            ans$par<-rep(NA,npar)
            ans$conv<-9999 # failed in run
            ans$gevals<-NA 
            ans$niter<-NA
         }
         if (ctrl$maximize) {
            ans$value= -ans$value
         }
      }  ## end if using newuoa
## --------------------------------------------
##      else if (meth == "DEoptim") {# Use DEoptim package
##        stop("This code for DEoptim has not been verified")
##        if (! is.null(mcontrol$maxit)) { 
##      mcontrol$maxfun<-mcontrol$maxit
##   } else {
##      mcontrol$maxfun<-500000 # ?? default at 091018, but should it be changed?!!
##   }
##        mcontrol$iprint<-mcontrol$trace
#        time <- system.time(ans <- try(DEoptim(fn=fn, lower=lower, upper=upper, control=mcontrol,...), silent=TRUE))[1]
##        if (class(ans)[1] != "try-error") {
##      ans$conv <- 0
##                if (ans$feval > mcontrol$maxfun) {
##         ans$conv <- 1 # too many evaluations
##                }
##           ans$fevals<-ans$feval
##           ans$gevals<-NA
##      ans$value<-ans$fval 
##        } else {
##      if (ctrl$trace>0) cat("bobyqa failed for current problem \n")
##             ans$conv<-9999
##           ans$fevals<-NA
##           ans$gevals<-NA
##      ans$value<-NA
##        }
##        ans$niter<-NULL
##      }  ## end if using DEoptim
## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      times[i] <- times[i] + time # Accumulate time for a single method (in case called multiple times)
      msg<-paste(meth," failed!")
      if (ans$conv == 9999) stop(msg)
      if (ctrl$trace>0) { cat("Post processing for method ",meth,"\n") }
      pars<-NULL # Unless length(pars)>0 we don't save anything later
      gradOK<-FALSE
      hessOK<-FALSE # to ensure these are set when kkt FALSE
      if ( ctrl$save.failures || (ans$conv < 1) ){  # Save the solution if converged or directed to save
         j <- j + 1 ## increment the counter for (successful) method/start case
         if (ctrl$trace && ans$conv==0) cat("Successful convergence! \n")  ## inform user we've succeeded
         # Testing final solution. Use numDeriv to get gradient and Hessian; compute Hessian eigenvalues
         gH<-NULL
         if ((ctrl$kkt || hessian) && (ans$conv != 9999)) {
            if (ctrl$trace>0) cat("Compute gradient and Hessian approximations at finish of ",method[i],"\n")
            gH<-gHgenb(ans$par, fn, gr=gr, hess=hess, bdmsk=bdmsk, lower=lower,
                 upper=upper, control=list(ktrace=ctrl$trace), ...) 
            # Note: Here we return to original functions supplied by user.
            # ??? 110703 maximize NOT accounted for.
            gradOK<-gH$gradOK
            hessOK<-gH$hessOK
            ngatend<-gH$ng
            nhatend<-gH$nH
         } # end test if hessian computed; note gradient is also computed
         ans$kkt1<-NA
         ans$kkt2<-NA
         if (gradOK) { # have gradient
            if (hessOK) { # have hessian
               akkt<- kktc(ans$par, ans$value, ngatend, nhatend, gH$nbm, 
			control=list(ktrace=ctrl$trace))
               ans$kkt1<-akkt$kkt1
               ans$kkt2<-akkt$kkt2
               if (ctrl$trace) {
                  cat("KKT results: gmax=",akkt$gmax,"  evratio=",akkt$evratio,
                          "  KKT1 & 2: ",akkt$kkt1,akkt$kkt2,"\n")
               }
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
            ans$systime <- time
            # Do we want more information saved?
            if (ctrl$trace>0) { 
               cat("Save results from method ",meth,"\n") 
            }
            ans.ret[[j]] <- ans  ## save the answer. [[]] indexes the CONTENTS of the list
            ans.ret[[j]]$method <- method[i] # and we tag on the method with the $ linker
            if (ctrl$trace>2) { cat("Assemble the answers\n") }
            #    attr(ans.ret, "CPU times (s)") <- times ## save the accumulated times 
            #    if (ctrl$trace>0) { cat("Done CPU times\n") }
            pars <- lapply(ans.ret, function(x) x$par)
            vals <- lapply(ans.ret, function(x) x$value)
            meths <- lapply(ans.ret, function(x) x$method)
            fevals<- lapply(ans.ret, function(x) x$fevals)
            gevals<- lapply(ans.ret, function(x) x$gevals)
            nitns <-  lapply(ans.ret, function(x) x$niter)
            convcode<- lapply(ans.ret, function(x) x$conv)
            kkt1<- lapply(ans.ret, function(x) x$kkt1)
            kkt2<- lapply(ans.ret, function(x) x$kkt2)
            xtimes<- lapply(ans.ret, function(x) x$systime)
         }  ## end post-processing of successful solution
         if (ctrl$follow.on) {
            par <- ans$par # save parameters for next method
            if (i < nmeth && ctrl$dowarn && ctrl$trace>0) cat("FOLLOW ON!\n")
         }
      } ## end loop over method (index i)
      if (length(pars) > 0) { # cannot save if no answers
         ansout <- data.frame(cbind(par=pars, fvalues=vals, method=meths, fns=fevals, grs=gevals, 
                        itns=nitns, conv=convcode, KKT1=kkt1, KKT2=kkt2, xtimes=xtimes))
         attr(ansout, "details") <- ans.ret
         # saveRDS(ansout,file="unsortedans")
         # save to file
         # sort by function value (DECREASING so best is last and 
         # follow.on gives natural ordering)
         # if (ctrl$trace>0) { cat("Sort results\n") }  # remove when working ?? 
         # print(ansout)
         if (ctrl$sort.result) { # sort by fvalues decreasing
            if(ctrl$maximize) ord <- order(as.numeric(ansout$fvalues))
            else ord <- rev(order(as.numeric(ansout$fvalues)))
            ansout <- ansout[ord, ]
         }
      } else {
      ansout<-NULL # no answer if no parameters
   } 
   ansout # return(ansout)
} ## end of optimx
