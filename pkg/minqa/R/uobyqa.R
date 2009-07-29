uobyqa.control <- function(rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000, wsize=NA)
  list(rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
uobyqa <- function(par, fn, control = uobyqa.control(), ...)
{
  n <- length(par) 
  if(n < 2)
    stop("uobyqa is not for optimization of single parameter.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- uobyqa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if (is.na(ctrl[["rhobeg"]]))
      ctrl[["rhobeg"]] <- max(1, max(abs(par))/2)   
  if (is.na(ctrl[["rhoend"]])) 
      ctrl[["rhoend"]] <- max(1.e-06, max(abs(par))/1e+06)
  w <- ( n**4 + 8*n**3 + 23*n**2 + 42*n + max(2*n**2 + 4, 18*n )) / 4
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")
 
  if (ctrl$maxfun < 10 * n^2) ctrl$maxfun = 10 * n^2
  out <- .Call("uobyqa_c", par, fn1, ctrl, new.env(), PACKAGE = "minqa")
  
  class(out) <- "uobyqa"
  out
}
