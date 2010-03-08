uobyqa.control <- function(rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000)
  list(rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun)
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

# rhobeg
  if(is.na(ctrl[["rhobeg"]])) {
      ctrl[["rhobeg"]] <- 0.2*max(abs(par))
      ctrl[["rhobeg"]]<-max(0.95, ctrl[["rhobeg"]]) # JN 090915??    
  }
# rhoend
  if(is.na(ctrl[["rhoend"]])) {
      ctrl[["rhoend"]]<-1.0e-6*ctrl[["rhobeg"]] # may want to change this.
  }

  w <- ( n**4 + 8*n**3 + 23*n**2 + 42*n + max(2*n**2 + 4, 18*n )) / 4
  ctrl[["wsize"]] <- w
  
  if (ctrl$maxfun < 10 * n^2) ctrl$maxfun = 10 * n^2
  out <- .Call("uobyqa_c", par, fn1, ctrl, new.env(), PACKAGE = "minqa")
  
  class(out) <- "uobyqa"
  out
}
