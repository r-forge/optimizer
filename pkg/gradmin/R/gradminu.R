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


  
if (control$trace > 2) {
    cat("control:")
    print(control)
}
  
npar <- length(x0)

# set up workspace
ctrl <- list(
  srchdirn = vmrf,
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
  status = "start",
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
  lsfail = FALSE, 
  lastsd = 0
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
#  w$lastsd <- 1 # for VM methods "last gradient"
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
  cat("method: ", deparse(substitute(w$srchdirn)),"\n")
  tmp <- readline(" ")
  if (w$srchdirn(w, msetup=TRUE,...) != 0) 
       stop(paste("Could not commence method ",deparse(substitute(w$srchdirn)),"\n"))
  
  repeat {
        if (w$trace > 0) {cat("Iteration ",w$niter,"(",w$ng,",",w$nf,")  Fval=",w$fbest,"\n")}
     if (terminate(w, ...)) break
     stp <- w$srchdirn(w,...) # no msetup 
     if (w$trace > 2) { 
         cat("srchdirn returned w$grd:")
         print(w$grd)
     }
     if (w$trace > 1) {
         cat("Search vector:")
         print(w$tdir)
     }
     if (w$trace > 2) {
        gprj <- as.numeric(crossprod(w$tdir, w$grd))
        cat("Gradient projection = ",gprj,"\n")
     }
     ## Do line search
     w$lsfail <- FALSE
     w$stsize<-w$lnsrch(w, ...)
     if (attr(w$stsize, "FAIL")) {
        if ((w$trace > 0) && (w$trace > 0)) cat("Line search fails\n")
        w$lsfail <- TRUE
     } else { # ls has succeeded in some way
       fval <- attr(w$stsize,"Fval")
      if (w$trace > 1) {cat(" step =", w$stsize,"  fval=",fval ," w$nf=",w$nf,"\n")}
      w$xn<-w$xb+w$stsize*w$tdir
      if (w$niter >= w$maxit) {
        print("Too many iterations. Failed to converge!")
        return(0)
      }
      w$xb <- w$xn
      w$fbest <- fval
      w$niter <- w$niter + 1
      if (w$watch) tmp <- readline("end iteration")   
     } # end else on lsfail
  } # end repeat
  out<-NULL # ensure cleared first, and then use structure above
# w$msg <- paste(deparse(substitute(w$srchdirn))," Apparently Successful", sep='')
  out <- list(w$xb, w$fbest, c(w$nf, w$ng, w$nh), w$convcode, w$msg, w$H)  #
  names(out) <- c("par", "value", "counts", "convergence", 
                  "message", "hessian")
  out
  }
