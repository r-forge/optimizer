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
  
  if(is.na(ctrl[["rhobeg"]]))
    ctrl[["rhobeg"]] <- abs(max(par) / 2)
  if(is.na(ctrl[["rhoend"]]))
    ctrl[["rhoend"]] <- abs(max(par) / 10e5)
  w <- ( n**4 + 8*n**3 + 23*n**2 + 42*n + max(2*n**2 + 4, 18*n )) / 4
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- (ctrl[["npt"]]+5)*(ctrl[["npt"]]+n)+3*n*(n+5)/2
  out <- .Call("uobyqa_c", unlist(par), fn1, ctrl, new.env(), PACKAGE = "minqa")
  
  class(out) <- "uobyqa"
  out
}
print.uobyqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("UOBYQA results\n")

    cat("parameter estimates:", toString(x$par), "\n")
    cat("function evaluations:", toString(x$feval), "\n")
    cat("objective function value:", toString(x$fval), "\n")
    
    invisible(x)
}
