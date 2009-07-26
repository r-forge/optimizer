newuoa.control <- function(npt = NA, rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000, wsize=NA)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun,
       wsize=wsize)
newuoa <- function(par, fn, control = newuoa.control(), ...)
{
  n <- length(par) 
  if(n < 2)
    stop("newuoa is not for optimization of single parameter.")
  fn1  <- function(par) fn(par, ...)
  ctrl <- newuoa.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }

  if(is.na(ctrl[["npt"]]))
    ctrl[["npt"]] <- length(par) * 2
  else if((ctrl[["npt"]] < n+2) || (ctrl[["npt"]] > (n+1)*(n+2)/2))
    stop("npt is not in [len(par)+2, (len(par)+1)*(len(par)+2)/2)] ") 
  if(is.na(ctrl[["rhobeg"]]))
    ctrl[["rhobeg"]] <- abs(max(par) / 2)
  if(is.na(ctrl[["rhoend"]]))
    ctrl[["rhoend"]] <- abs(max(par) / 10e5)
  w <- ( ctrl[["npt"]]+13)*( ctrl[["npt"]]+n)+3*n*(n+3)/2
  if(is.na(ctrl[["wsize"]]))
    ctrl[["wsize"]] <- w
  else if(ctrl[["wsize"]] < w) stop("wsize is not large enough.")
  
  out <- .Call("newuoa_c", unlist(par), fn1, ctrl, new.env(), PACKAGE = "minqa")
  
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
