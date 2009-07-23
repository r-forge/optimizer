newuoa.control <- function(npt = NA, rhobeg = 1e-2, rhoend = 1e-4,
                           iprint = 3, maxfun=10000, wsize=10000)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
newuoa <- function(par, fn, control = newuoa.control(), ...)
{
  if(length(par) < 2)
    warning("newuoa is not for optimization of single parameter.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- newuoa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(is.na(ctrl[["npt"]]))
    ctrl[["npt"]] <- length(par) + 2
  out <- .Call("newuoa_c", par, fn1, ctrl, new.env(), PACKAGE = "minqa")
  
  class(out) <- "newuoa"
  out
}
print.newuoa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("NEWUOA results\n")

    cat("parameter estimates:", toString(x$par), "\n")
    cat("function evaluations:", toString(x$feval), "\n")
    cat("objective function value:", toString(x$fval), "\n")
    
    invisible(x)
}
