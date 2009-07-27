bobyqa.control <- function(npt = NA, rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000, wsize=NA)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
bobyqa <- function(par, fn, xl, xu, control = bobyqa.control(), ...)
{

  n <- length(par) 
  if(n < 2)
    stop("bobyqa is not for optimization of single parameter.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- bobyqa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(is.na(ctrl[["npt"]]))
    ctrl[["npt"]] <- min(n * 2, n+2)
# Note: changed from n+1 by JN 090726
  else if((ctrl[["npt"]] < n+2) || (ctrl[["npt"]] > (n+1)*(n+2)/2))
    stop("npt is not in [len(par)+2, (len(par)+1)*(len(par)+2)/2)] ") 
  if(ctrl[["npt"]] > (2*n -1) )
    warning("Setting 'npt' larger than 2 * length(par)+1 not recommended.")
  if(is.na(ctrl[["rhobeg"]])) {
    if(all(is.finite(xu-xu))) 
      ctrl[["rhobeg"]] <- max(xu-xl)
    else 
      ctrl[["rhobeg"]] <-  if(all(abs(par)<2)) 1 else max(abs(par) / 2)
  }
  if(is.na(ctrl[["rhoend"]])) {
    if(all(is.finite(xu-xu))) 
      ctrl[["rhoend"]] <- max(xu-xl) / 10e5
    else 
      ctrl[["rhoend"]] <- max(abs(par) / 10e5) 
  }
  if(ctrl[["rhobeg"]] < ctrl[["rhoend"]] ||
     any(c(ctrl[["rhobeg"]],ctrl[["rhoend"]]) < 0))
    ## may not need this (Inf and -Inf seem to work as bound defaults)
    ##||
     ##any(xu-xl > 2 * ctrl[["rhobeg"]]))
    stop("rhobeg and rhoend must be positive with rhoend no greater than
       rhobeg.") 
  if(all(is.finite(xu-xu)) && any(xu-xl > 2 * ctrl[["rhobeg"]]))
    stop("All of the differences xu-xl must be less than 2*rhobeg.") 

  w <- (ctrl[["npt"]]+5)*(ctrl[["npt"]]+n)+3*n*(n+5)/2
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")
  
  out <- .Call("bobyqa_c", par, xl, xu, fn1, ctrl, new.env(),
               PACKAGE = "minqa")
  
  class(out) <- "bobyqa"
  out
}
print.bobyqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("BOBYQA results\n")

    cat("parameter estimates:", toString(x$par), "\n")
    cat("function evaluations:", toString(x$feval), "\n")
    cat("objective function value:", toString(x$fval), "\n")
    
    invisible(x)
}
