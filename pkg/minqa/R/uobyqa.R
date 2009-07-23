uobyqa.control <- function(rhobeg = 1e-2, rhoend = 1e-4,
                           iprint = 3, maxfun=10000, wsize=10000)
  list(rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
uobyqa <- function(par, fn, control = uobyqa.control(), ...)
{
  if(length(par) < 2)
    warning("uobyqa is not for optimization of single parameter.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- uobyqa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
 
  out <- .Call("uobyqa_c", par, fn1, ctrl, new.env(), PACKAGE = "minqa")
  
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
