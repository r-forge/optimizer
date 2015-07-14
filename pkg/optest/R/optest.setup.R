optest.setup <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

#  cat("optest.setup, control is \n")
#  print(control)
#  tmp <- readline("continue")

# To return in optcfg: fname, npar ??, method, ufn, ugr, ctrl

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
#   if (is.null(control)) {
   ctrl <- ctrldefault(npar)
#   } else { ctrl <- control }
# substitute the appropriate information according to user spec.
   if (! is.null(control)) { 
      ncontrol <- names(control)
      nctrl <- names(ctrl)
      for (onename in ncontrol) {
         if (onename %in% nctrl) {
            ctrl[onename] <- control[onename]
         } 
      }
   } else { 
      cat("setting control defaults \n")
      ctrl <- ctrldefault(npar)
   }

## Now we need to define user functions that are w r t SCALED parameters

# See if we should do kkt check (default is true)
# control$kkt defined, then use that value
  if (is.null(control$kkt)) {  # Only adjust kkt control if user has not specified it
    if (ctrl$kkt) { # if kkt is on, turn it off for large matrices (default 150504 is TRUE)
      if (is.null(gr) || is.character(gr)) { # no analytic gradient
         if (npar >  ctrl$grcheckfnog) {
           ctrl$kkt=FALSE # too much work when large number of parameters
           if (ctrl$trace>0) cat("gr approx., npar > ",ctrl$grcheckfnog,", ctrl$kkt set FALSE\n")
         }
      } else {
         if (npar > ctrl$grcheckfwithg) {
            ctrl$kkt=FALSE # too much work when npar large, even with analytic gradient
            if (ctrl$trace>0) cat("Have gr, npar > ",ctrl$grcheckfwithg,", ctrl$kkt set FALSE\n")
         }
      }
    }
  }

# Set the scaled parameters
   optcfg$spar <- par / ctrl$parscale # we DIVIDE by the scalings
   optcfg$lower <- lower / ctrl$parscale # bounds adjust
   optcfg$upper <- upper / ctrl$parscale

# reset the function if we are maximizing / using a gradient approx.
  ctrl$maximizeorig <- ctrl$maximize
  ctrl$maximize <- FALSE # We always minimize internally
  if (ctrl$usenumDeriv) { 
     gr <- "grnd" 
  } else if (is.null(gr)) { gr <- ctrl$defgrapprox } # use default gradient approximation
  defscale <- rep(1, npar)
  if (ctrl$maximizeorig) {
     ufn <- function(spar, userfn=fn, pscale=defscale, ...){
# note pscale is a vector
        result <- (-1) * fn(spar*pscale, ...)
     }
     if ( is.character(gr) ) {
       # Convert string to function call, assuming it is a numerical gradient function
       ugr <- function(par=par, userfn=fn, pscale=defscale, ...){
           result <- do.call(gr, list(par, userfn, ...))
           result <- (-1)*result*pscale
       }
       uhess <- NULL # Do NOT define hessian when approximating gradient
     } else { ugr<- function(par=par, userfn=fn, pscale=defscale, ...) {
                 result <- (-1)*gr(par,...)*pscale
              }
              if (! is.null(hess) ) {
                 uhess <- function(par=par, userfn=fn, pscale=defscale, ...) {
                    result <- (-1) * (diag(pscale) %*% hess(par*pscale, ...) %*% diag(pscale))
                 }
              }   
     }
  } else { # not maximizing
       ufn <- function(spar, userfn=fn, pscale=defscale, ...){
           val <- fn(spar * pscale, ...)
       }
       if ( is.character(gr) ) {
         # Convert string to function call, assuming it is a numerical gradient function
         ugr<-function(par=par, userfn=fn, pscale=defscale, ...){
           result <- do.call(gr, list(par*pscale, userfn, ...))*pscale
         }
         uhess <- NULL # Do NOT define hessian when approximating gradient
       } else { 
           ugr <- function(par=par, userfn=fn, pscale=defscale, ...){ # to satisfy formal arguments
               result <- gr(par*pscale, ...)*pscale
#              result <- do.call(gr, list(par,...))*pscale
            }
           if (is.null(hess)) {
              uhess <- NULL
           } else {
              uhess <- function(par=par, userfn=fn, pscale=defscale, ...) {
                  (diag(pscale) %*% hess(par*pscale, ...) %*% diag(pscale))
              }
           }
       }
  } # set functions if MAXIMIZEing
  optcfg$usenumDeriv <- FALSE # JN130703 ?? Because we have already set ugr
  optcfg$ufn <- ufn
  optcfg$ugr <- ugr
  optcfg$uhess <- uhess

  optcfg$ctrl <- ctrl
  optcfg$method <- method
  optcfg # return the structure
   openv <<- list2env(optcfg, parent=emptyenv())
} ## end of optimx.setup
