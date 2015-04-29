optimx.setup <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

# To return in optcfg: fname, npar ??, method, ufn, ugr, ctrl, have.bounds

# Get real name of function to be minimized
  fname<-deparse(substitute(fn))
  if (!is.null(control$trace) && control$trace>0) {
	cat("fn is ",fname,"\n")
  }
  optcfg <- list()
  optcfg$fname <- fname

# Only one ref to parameters -- to get npar here
  npar <- length(par) # !! NOT CHECKED in case par not well-defined
  optcfg$npar <- npar

##??140828 need default and alternate sets of controls and methods

## We have npar, so can compute various controls like maxfeval

# Set control defaults # ?? does not include some of the tests and maxit etc.
    ## for different methods

    ctrl.default <- list(
        acctol = 0.0001, 
	all.methods=FALSE,
	badval=(0.5)*.Machine$double.xmax,
	dowarn=TRUE, 
        eps = 1e-07, 
	follow.on=FALSE, 
        grcheckfwithg=500,
        grcheckfnog=50,
	kkt=TRUE,
	kkttol=0.001,
	kkt2tol=1.0E-6,
	maximize=FALSE,
        maxit=500*round(sqrt(npar+1)),
	maxfeval=5000*round(sqrt(npar+1)),
        reltest=100.0,
	save.failures=TRUE,
	scaletol=3, 
	starttests=TRUE,
        stepredn=0.2,
        stopbadupdate=FALSE,
	trace=0,
        usenumDeriv=FALSE
    ) 
    
# Control set for optimx version 2013.8.6
# Note: same as default -- must be in maxit etc. ??!!
   ctrl.2013 <- list(
        acctol = 0.0001, 
	all.methods=FALSE,
	badval=(0.5)*.Machine$double.xmax,
	dowarn=TRUE, 
        eps = 1e-07, 
	follow.on=FALSE, 
        grcheckfwithg=500,
        grcheckfnog=50,
	kkt=TRUE,
	kkttol=0.001,
	kkt2tol=1.0E-6,
	maximize=FALSE,
        maxit=500,
	maxfeval=5000*round(sqrt(npar+1)),
        reltest=100.0,
	save.failures=TRUE,
	scaletol=3, 
	starttests=TRUE,
        stepredn=0.2,
        stopbadupdate=FALSE,
	trace=0,
        usenumDeriv=FALSE
   )

   if (is.null(control)) { ctrl <- ctrl.default }
   else { if (! is.null(control$ctrlset) ) { # set to another set of controls
            if (control$ctrlset=="default") {ctrl <- ctrl.default}
            else if (control$ctrlset=="2013") {ctrl <- ctrl.2013}
                 else { stop("BAD END: control set not defined") }
        } else { ctrl <- ctrl.default }
   }
   

# Note that we do NOT want to check on the names, because we may introduce 
#    new names in the control lists of added methods
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
# However, we do want to substitute the appropriate information. 
# removed copy of hessian to control$kkt
##?? Clean this up??
    ncontrol <- names(control)
    nctrl <- names(ctrl)
#    cat("names in control:")
#    print(ncontrol)
#    cat("names in ctrl:")
#    print(nctrl)
    for (onename in ncontrol) {
       if (onename %in% nctrl) {
           ctrl[onename] <- control[onename]
       } else {
           ctrl[onename] <- control[onename]
       }
    }
## ?? parametrize the 50 and 500
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
        if (is.null(gr)) {
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
    cat("optimx.setup ctrl out:")
    print(ctrl)


    optcfg$ctrl <- ctrl
# reset the function if we are maximizing
  ufn <- fn
  ugr <- gr
  uhess <- hess
  if ((! is.null(control$maximize)) && control$maximize ) { 
        cat("Maximizing -- use negfn and neggr\n")
        if (! is.null(control$fnscale)) { 
 		stop("Mixing controls maximize and fnscale is dangerous. Please correct.")
        } # moved up 091216
        optcfg$ctrl$maximize<-TRUE
        ufn <- function (par, ...) { # negate the function for maximizing
	   val<-(-1.)*fn(par,...)
        } # end of ufn = negfn
        if (! is.null(gr)) { 
           ugr <- function(par, userfn=ufn, ...) {
               gg <- (-1)*gr(par, ...)
           }
        } else { ugr <- NULL } # ensure it is defined
        if (! is.null(hess) ) {
           uhess <- function(par, ...) {
               hh <- (-1)*hess(par, ...)
           }
        } else { uhess <- NULL } # ensure it is defined
  } else { 
     optcfg$ctrl$maximize <- FALSE # ensure defined
  } # define maximize if NULL
  optcfg$usenumDeriv <- FALSE # JN130703
  if (is.null(gr) && ctrl$usenumDeriv) {
     if (ctrl$dowarn) warning("Replacing NULL gr with 'numDeriv' approximation")
     optcfg$usenumDeriv<-TRUE
     ugr <- function(par, userfn=ufn, ...) { # using grad from numDeriv
        tryg<-grad(userfn, par, ...)
     } # Already have negation in ufn if maximizing
  }
  if (is.character(gr)) {
     if (ctrl$dowarn) warning("Replacing NULL gr with '",gr,"' approximation")
     ugr <- function(par, userfn=ufn, ...) { # using grad from numDeriv
        tryg <- do.call(gr, list(par, userfn, ...))
     } # Already have negation in ufn if maximizing
  }
  optcfg$ufn <- ufn
  optcfg$ugr <- ugr
  optcfg$uhess <- uhess

# Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) { have.bounds<-TRUE # set this for convenience
  } else { have.bounds <- FALSE }
  optcfg$have.bounds <- have.bounds

# List of methods in base or stats, namely those in optim(), nlm(), nlminb()
  basemeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb")
# SANN has no termination for optimality, only a maxit count for
#   the maximum number of function evaluations; remove DEoptim for now -- not useful 
#   for smooth functions. Code left in for those who may need it.
# List of methods in this packages
   oxmeth <- c("lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin")
   pkgmeth <- c("spg", "ucminf", "newuoa", "bobyqa", "nmkb", "hjkb")

# Now make sure methods loaded
   allmeth <- c(basemeth, oxmeth, pkgmeth) 

   goodpkg <- (requireNamespace("BB", quietly=TRUE) && 
             requireNamespace("ucminf", quietly=TRUE) && 
             requireNamespace("minqa", quietly=TRUE) && 
             requireNamespace("dfoptim", quietly=TRUE) && 
             requireNamespace("setRNG", quietly=TRUE)) 

   if (! goodpkg) stop("Installation missing one of BB, ucminf, minqa, dfoptim, setRNG")
  
   bdsmeth<-c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", 
     "Rtnmin", "bobyqa", "nmkb", "hjkb", "lbfgsb3")
  # Restrict list of methods if we have bounds
  if (any(is.finite(c(lower, upper)))) allmeth <- bdsmeth
  if (("All" %in% method) || ("ALL" %in% method)) stop("To specify all methods, use 'all' (lower case)")
  if (("all" %in% method) && (length(method) > 1)) stop("If method='all', method can have ONLY that element")
  if (method == "all") ctrl$all.methods=TRUE # use ctrl$all.methods as the main mechanism
  if (! ctrl$all.methods) { # Changes method vector!
      if ((length(method) == 1) && (method=="all")) ctrl$all.methods <- TRUE
  }             
  if (ctrl$all.methods) { # Changes method vector!
	method <- allmeth
        if (ctrl$trace > 0) {
		cat("all.methods is TRUE -- Using all available methods\n")
		print(method)
	}
  } 

  # Partial matching of method string allowed
  # avoid duplicates here
  # 2011-1-17 JN: to set L-BFGS-B
  if (is.null(method)) {
     if (is.null(gr)) { method <- "nmkb" }
     else { 
        if (have.bounds) { method <- "L-BFGS-B" }
        else {method <- "BFGS" }
     }
  }
  method <- try(unique(match.arg(method, allmeth, several.ok=TRUE) ),silent=TRUE)
  if (class(method)=="try-error") stop("optimx: No match to available methods")
  nmeth <- length(method)
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
  optcfg$method <- method
  optcfg # return the structure
} ## end of optimx.setup
