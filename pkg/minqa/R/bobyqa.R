bobyqa.control <- function(npt = NA, rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000, wsize=NA,
                           obstop=TRUE, force.start=FALSE)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize, obstop=obstop, force.start=force.start)
########################################################################
bobyqa <- function(par, fn, lower=-Inf, upper=Inf, control = bobyqa.control(), ...)
{
  n <- length(par) 
  if(n < 2)
    stop("bobyqa is not for optimization of single parameter.")
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  if(length(lower)!=n || length(upper)!=n)
      stop("Bounds are the wrong length.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- bobyqa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(is.na(ctrl[["npt"]]))
    ctrl[["npt"]] <- min(n+6, 2*n+1) # Follow Powell's lead
# Note: changed from n+1 by JN 090726
  else if((ctrl[["npt"]] < n+2) || (ctrl[["npt"]] > (n+1)*(n+2)/2))
    stop("npt is not in [len(par)+2, (len(par)+1)*(len(par)+2)/2)] ") 
  if(ctrl[["npt"]] > (2*n + 1) )
    warning("Setting 'npt' larger than 2 * length(par)+1 not recommended.")
  if (ctrl[["iprint"]]>1) cat("npt is set at ",ctrl[["npt"]]," and n=",n,"\n")
  if(is.na(ctrl[["rhobeg"]])) {
    if(all(is.finite(upper-lower))) {
      if ( any(upper <= lower) ) stop("Overlapping bounds")
      ctrl[["rhobeg"]] <- 0.2*mean(upper-lower)
    }
    else 
      ctrl[["rhobeg"]] <- 0.2*max(abs(par))
  }
  if(is.na(ctrl[["rhoend"]])) {
      ctrl[["rhoend"]]<-1.0e-6*ctrl[["rhobeg"]] # may want to change this.
  }
  if (ctrl[["iprint"]]>1) cat("RHOBEG =",ctrl[["rhobeg"]],"  RHOEND=",ctrl[["rhoend"]],"\n")
  if(ctrl[["rhobeg"]] < ctrl[["rhoend"]] ||
     any(c(ctrl[["rhobeg"]],ctrl[["rhoend"]]) < 0))
    ## may not need this (Inf and -Inf seem to work as bound defaults)
    ##||
     ##any(upper-lower > 2 * ctrl[["rhobeg"]]))
    stop("rhobeg and rhoend must be positive with rhoend no greater than rhobeg.") 
  if(all(is.finite(upper-lower)) && any(upper-lower < 2 * ctrl[["rhobeg"]])) {
    stop("All of the differences upper-lower must be >= 2*rhobeg.") 
  }
  w <- (ctrl[["npt"]]+5)*(ctrl[["npt"]]+n)+3*n*(n+5)/2
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")

  # Out of bounds check.
  if (all(is.finite(lower)) ) {
      if (any(par < lower) ) {
          msg<-"Some parameters out of bounds LOW"
          if (ctrl[["obstop"]]) {
             stop(msg)
          } else {
             warning(msg)
             par[which(par < lower)] <- lower
             warning("Some parameters adjusted to the nearest bound")
	  }
      }          
  }
  if (all(is.finite(upper)) ) {
      if (any(par > upper) ) {
          msg<-"Some parameters out of bounds HIGH"
          if (ctrl[["obstop"]]) {
             stop(msg)
          } else {
             warning(msg)
             par[which(par > upper)] <- upper
             warning("Some parameters adjusted to the nearest bound")
	  }
      }          
  }

  if (all(is.finite(upper)) && all(is.finite(lower)) && all(par >= lower) && all(par <= upper) ) {
     # everything in bounds. Decide if we move parameters.
     if (ctrl[["iprint"]]>1) cat("ctrl[[force.start]] = ", ctrl[["force.start"]],"\n")
     if (! ctrl[["force.start"]] ) {
        if (any( (par - lower) < ctrl[["rhobeg"]]) ){
           par[which( (par - lower) < ctrl[["rhobeg"]] )] <- lower + ctrl[["rhobeg"]]
           warning("Some parameters adjusted away from lower bound")
        }
        if (any( (upper - par) < ctrl[["rhobeg"]]) ) {
           par[which( (upper - par) < ctrl[["rhobeg"]] )] <- upper - ctrl[["rhobeg"]]
           warning("Some parameters adjusted away from upper bound")
        }
     }
  }

  if (ctrl$maxfun < 10 * n^2)
     warning("'maxfun' less than 10*length(par)^2 not recommended.")

  out <- .Call("bobyqa_c", unlist(par), lower, upper, fn1, ctrl, new.env(),
               PACKAGE = "minqa")
  
  class(out) <- "bobyqa"
  out
}
print.bobyqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("bobyqa results\n")

    cat("parameter estimates:", toString(x$par), "\n")
    cat("function evaluations:", toString(x$feval), "\n")
    cat("objective function value:", toString(x$fval), "\n")
    
    invisible(x)
}
