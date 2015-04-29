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
  optcfg$npar <- npar ##?? Need to put this stuff in environment??

##??140828 need default and alternate sets of controls and methods

## We have npar, so can compute various controls like maxfeval

# Set control defaults # ?? does not include some of the tests and maxit etc.
    ## for different methods
   ctrl <- ctrldefault(npar)

# substitute the appropriate information according to user spec.

   if (! is.null(control)) { 
      ncontrol <- names(control)
      nctrl <- names(ctrl)
      for (onename in ncontrol) {
         if (onename %in% nctrl) {
            ctrl[onename] <- control[onename]
         } 
      }
   }

# Check if we should do gradient check
    if (is.null(control$kkt)) { # turn off kkt for large matrices
      ctrl$kkt<-TRUE # default it to compute KKT tests
      if (is.null(gr)) { # no analytic gradient
         if (npar >  ctrl$grcheckfnog) {
           ctrl$kkt=FALSE # too much work when large number of parameters
           if (ctrl$trace>0) cat("gr NULL, npar > ",ctrl$grcheckfnog,", ctrl$kkt set FALSE\n")
         }
      } else {
         if (npar > ctrl$grcheckfwithg) {
            ctrl$kkt=FALSE # too much work when npar large, even with analytic gradient
            if (ctrl$trace>0) cat("gr NULL, npar > ",ctrl$grcheckfwithg,", ctrl$kkt set FALSE\n")
         }
      }
    } else { # kkt is set
      if (control$kkt) {
        if (is.null(gr)) {
           if (npar > ctrl$grcheckfnog) {
             if ((ctrl$trace>0) && ctrl$dowarn) 
                warning("Computing hessian for gr NULL, npar > ", ctrl$grcheckfnog,", can be slow\n")
           }
        } else {
           if (npar > ctrl$grcheckfwithg) {
             if ((ctrl$trace>0) && ctrl$dowarn) 
               warning("Computing hessian with gr code, npar > ",ctrl$grcheckfwithg,", can be slow\n")
           }
        }
      }
    }
##    cat("optimx.setup ctrl out:")
##    print(ctrl)

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
