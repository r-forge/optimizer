bobyqa.control <- function(npt = NA, rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000, wsize=NA)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
bobyqa <- function(par, fn, lower=-Inf,
                   upper=Inf, control = bobyqa.control(), ...)
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
  xl <- lower
  xu <- upper
  if (length(xl) == 1) xl <- rep(xl, n)
  if (length(xu) == 1) xu <- rep(xu, n)
  if(length(xl)!=n || length(xu)!=n)
    stop("Bounds are the wrong length.")
  if(is.na(ctrl[["npt"]]))
    ctrl[["npt"]] <- n+2
  else if((ctrl[["npt"]] < n+2) || (ctrl[["npt"]] > (n+1)*(n+2)/2))
    stop("npt is not in [len(par)+2, (len(par)+1)*(len(par)+2)/2)] ") 
  if(ctrl[["npt"]] > (2*n +1) )
    warning("Setting 'npt' larger than 2 * length(par)+1 not recommended.")
  if (is.na(ctrl[["rhobeg"]])) {
        if (all(is.finite(xu - xl))) 
          ctrl[["rhobeg"]] <- min(xu - xl)/2
        else ctrl[["rhobeg"]] <- max(1, max(abs(par))/2)
      }
  if (is.na(ctrl[["rhoend"]])) {
    if (all(is.finite(xu - xl))) 
      ctrl[["rhoend"]] <- max(xu - xl)/1e+06
    else ctrl[["rhoend"]] <- max(1.e-06, max(abs(par))/1e+06)
  }
  if(ctrl[["rhobeg"]] < ctrl[["rhoend"]] ||
     any(c(ctrl[["rhobeg"]],ctrl[["rhoend"]]) < 0))
    ## may not need this (Inf and -Inf seem to work as bound defaults)
    ##||
     ##any(xu-xl > 2 * ctrl[["rhobeg"]]))
    stop("rhobeg and rhoend must be positive with rhoend no greater than
       rhobeg.") 
  if(all(is.finite(xu-xl)) && any(xu-xl < 2 * ctrl[["rhobeg"]]))
    stop("None of the differences upper-lower can be less than 2*rhobeg.")

  w <- (ctrl[["npt"]]+5)*(ctrl[["npt"]]+n)+3*n*(n+5)/2
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")
  if (ctrl$maxfun < 10 * n^2) ctrl$maxfun = 10 * n^2
  out <- .Call("bobyqa_c", par, xl, xu, fn1, ctrl, new.env(),
               PACKAGE = "minqa")
  
  class(out) <- "bobyqa"
  out
}
