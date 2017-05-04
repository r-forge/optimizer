gradminu<-function(x0, fn, gr, hess, lower = NULL, upper = NULL, 
      control=list(),...) {
## General (unconstrained) gradient minimizer 
##
##Input
##       - fn is the function we wish to minimize
##?? fixup documentation here??
##       - x0 is the initial value
##       - ... is data used in the function fn
##Output (list) -- matching optimr() output
#  par	  The best set of parameters found.
#  value	The value of ‘fn’ corresponding to ‘par’.
#  counts	A 3-element integer vector giving the number of calls to 
#      ‘fn’ and ‘gr’ and 'hess' respectively. 
#  convergence	An integer code. ‘0’ indicates successful completion
#  message	 A character string giving any additional information 
#  returned by the optimizer, or ‘NULL’.
#  hessian	Either NULL or the Hessian matrix at the parameters par.

# Need a way to pass any object defined at this level into routines called
# But not necessarily to define them WITHIN this routine

## ?? How to put in name of solver, lnsrch for display
## ?? Need a string/list to specify method with parameters, and output with result

terminate <- function(w, ...){
  if (w$trace > 2) cat("Termination test:")    
  halt <- FALSE # default is keep going
  # tests on too many counts??
  if (niter >= w$maxit) {
    if (w$trace > 0) cat("Too many (",niter," iterations\n")
    halt <- TRUE
    w$convcode <- 1
    break    
  }
  if (w$nf >= w$maxfevals) {
    w$msg <- paste("Too many (",w$nf,") function evaluations")
    if (w$trace > 0) cat(w$msg,"\n")
    halt <- TRUE
    w$convcode <- 1 # ?? value
    break
  }
  #    if (w$ng > w$maxgevals){} # not implemented
  #    if (w$nh > w$maxhevals){} # not implemented
  gmax <- max(abs(w$grd))
  if (gmax <= w$epstol*w$f0) {
    w$msg <- paste("Small gradient norm",gmax)
    if (w$trace > 0) cat(w$msg,"\n")
    halt <- TRUE
    w$convcode <- 0 # OK
  }
  halt
}  
  
if (control$trace > 2) {
    cat("control:")
    print(control)
}
  
npar <- length(x0)

# set up workspace
ctrl <- list(
  minmeth = vmrf,
  lnsrch = lsback,
  solver = solve,
  trace = 0,
  maxit = 500,
  maxfevals = npar*500,
  laminc = 10,
  lamdec = 0.4,
  lamstart = 0.01,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  stepdec = 0.2, 
  stepmax = 5,
  stepmin = 0,
  smult = 1.5,
  fmult = -0.25,
  offset = 100.0,
  defstep=1,
  bigval = .Machine$double.xmax*0.01,
  nf = 0,
  ng = 0,
  nh = 0,
  niter = 0,
  npar = npar,
  watch = FALSE,
  f0 = NA,
  stopiter = FALSE,
  convcode = NA,
  msg = "??",
  lastng = 0
)  

ncontrol <- names(control)
nctrl <- names(ctrl)
for (onename in ncontrol) {
  if (onename %in% nctrl) {
    ctrl[onename]<-control[onename]
  }
}

w <- list2env(ctrl) # Workspace
w$fn <- fn
w$gr <- gr
w$hess <- hess


  w$lambda<-w$lamstart ## ?? do better
  w$lastng <- 1 # for VM methods "last gradient"
  w$niter <- 1
  w$xb <- x0 # best so far
  w$fbest <- fn(w$xb, ...)
  w$nf <- w$nf + 1
  if ((class(w$fbest) == "try-error") | is.na(w$fbest) | 
                 is.null(w$fbest) | is.infinite(w$fbest)) {
    w$msg <- "Initial point gives inadmissible function value"
    w$convcode <- 20 # ?? choice of settings
    if (trace > 0)  cat(msg, "\n") 
    out <- list(w$xb, w$bigval, c(w$nf, 0, 0), w$convcode, w$msg)  #
    names(out) <- c("par", "value", "counts", "convergence", 
                    "message")
# ?? add Hess etc??
    return(out)
  }
  
  if (w$trace > 1) cat("Initial Fval =",w$fbest,"\n")
  w$f0 <- w$fbest # for scaling
  # Initialize
  if (w$minmeth(w, msetup=TRUE) != 0) 
       stop(paste("Could not commence method ",deparse(substitute(w$minmeth)),"\n"))
  cat("returned w$grd:")
  print(w$grd)
  repeat {
        if (w$trace > 0) {cat("Iteration ",w$niter,"  Fval=",w$fbest,"\n")}
     stp <- w$minmeth(w,...) # no msetup
     if (w$trace > 1) {
         cat("Search vector:")
         print(stp)
     }
     if (w$trace > 1) {
        gprj <- as.numeric(crossprod(stp, w$grd))
        cat("Gradient projection = ",gprj)
     }
     ## Do line search
     w$gvl<-w$lnsrch(w, ...)
     if (attr(w$gvl, "FAIL")) {
        if ((w$trace > 0) && (w$trace > 0)) cat("Line search fails\n")
        break
     }
     fval <- attr(w$gvl,"Fval")
##    cat("after lnsrch w$nf now ",w$nf,"\n")
    if (w$trace > 2) {cat(" step =", w$gvl,"  fval=",fval ," w$nf=",w$nf,"\n")}
    xn<-xb+w$gvl*stp
    if (niter >= w$maxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
    xb <- xn
    fbest <- fval
    niter <- niter + 1
    if (w$watch) tmp <- readline("end iteration")   
  } # end repeat
  out<-NULL # ensure cleared first, and then use structure above
  w$msg <- "Apparently Successful"
  out <- list(w$xb, w$fbest, c(w$nf, w$ng, w$nh), w$convcode, w$msg, w$H)  #
  names(out) <- c("par", "value", "counts", "convergence", 
                  "message", "hessian")
  out
  }
