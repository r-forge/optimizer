newuoa.control <- function(npt = NA, rhobeg = NA, rhoend = NA,
                           iprint = 0, maxfun=10000)
  list(npt = npt, rhobeg=rhobeg, rhoend=rhoend, iprint=iprint, maxfun=maxfun)
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
    ctrl[["npt"]] <- n+2
# Note: changed from n+1 by JN 090726
  else if((ctrl[["npt"]] < n+2) || (ctrl[["npt"]] > (n+1)*(n+2)/2))
    stop("npt is not in [len(par)+2, (len(par)+1)*(len(par)+2)/2)] ") 
#  if(is.na(ctrl[["rhobeg"]]))
#    ctrl[["rhobeg"]] <-  max(1, max(abs(par))/2)
#  if(is.na(ctrl[["rhoend"]]))
#    ctrl[["rhoend"]] <-  max(1.e-06, max(abs(par))/1e+06)
# rhobeg
  if(is.na(ctrl[["rhobeg"]])) {
      ctrl[["rhobeg"]] <- 0.2*max(abs(par))
      ctrl[["rhobeg"]]<-max(0.95, ctrl[["rhobeg"]]) # JN 090915??    
  }
# rhoend
  if(is.na(ctrl[["rhoend"]])) {
      ctrl[["rhoend"]]<-1.0e-6*ctrl[["rhobeg"]] # may want to change this.
  }


  w <- ( ctrl[["npt"]]+13)*( ctrl[["npt"]]+n)+3*n*(n+3)/2
  ctrl[["wsize"]] <- w
  
  if (ctrl$maxfun < 10 * n^2) ctrl$maxfun = 10 * n^2
  out <- .Call("newuoa_c", par, fn1, ctrl, new.env(), PACKAGE = "minqa")
  
  class(out) <- "newuoa"
  out
}
