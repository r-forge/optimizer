gradminu<-function(x0, fn, gr, hess, lower = NULL, upper = NULL, 
      control=list(),...) {
## General (unconstrained) gradient minimizer 
##
##Input
##       - fn is the function we wish to minimize
##?? fixup documentation here??
##       - x0 is the initial value
##       - ... is data used in the function fn
##Output (list) -- need to match optim() output!! ???
##       - xs is the value at the minimum
##       - fv is the fn evaluated at xs
##       - grd is the gradient
##       - Hess is the Hessian
##       - niter is the number of interations needed (gradient and Hessian evals).
##       - add fevals??, other reports

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
  watch = FALSE
)  

ncontrol <- names(control)
nctrl <- names(ctrl)
for (onename in ncontrol) {
  if (onename %in% nctrl) {
    ctrl[onename]<-control[onename]
  }
}

w <- list2env(ctrl) # Workspace

  lambda<-w$lamstart ## ?? do better
  niter <- 1
  xb <- x0 # best so far
  fbest <- fn(xb, ...)
  w$nf <- w$nf + 1
  if (w$trace > 0) cat("Initial Fval =",fbest,"\n")
  f0 <- fbest # for scaling
#  cat("w$nf now ",w$nf,"\n")
##  getH <- TRUE
  repeat {
##    if (getH) {
        if (w$trace > 0) {cat("Iteration ",niter,"  Fval=",fbest,"\n")}
        grd<-gr(xb,...)
        w$ng <- w$ng + 1
        H<-hess(xb,...)
        w$nh <- w$nh + 1
##    }
##    cat("Termination test:")    
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter >= w$maxit) {
        if (w$trace > 0) cat("Too many (",niter," iterations\n")
        halt <- TRUE
        convcode <- 1
        break
    }
    cat(" w$nf=",w$nf,"  w$maxfevals=",w$maxfevals,"\n")
    if (w$nf >= w$maxfevals) {
      cat("Stopping\n")
      if (w$trace > 0) cat("Too many ",w$nf," function evaluations\n")
      halt <- TRUE
      convcode <- 91 # ?? value
      break
    }
    #    if (w$ng > w$maxgevals){} # not implemented
    #    if (w$nh > w$maxhevals){} # not implemented
    gmax <- max(abs(grd))
    if (gmax <= w$epstol*f0) {
      if (w$trace > 0) cat("Small gradient norm",gmax,"\n")
      halt <- TRUE
      convcode <- 0 # OK
      break
    }
    stp<-try(w$solver(H, -grd))
    if (class(stp) == "class-error") {
          stop("Failure of default solve of Newton equations")
    }
##    } else if (w$solver == "marquardt") {
##       Haug<-H + (diag(H)+1.0)*lambda # To avoid singularity
##       stp <- solve(Haug, -grd)
    
    if (w$trace > 0) {
         cat("Search vector:")
         print(stp)
    }

    gprj <- as.numeric(crossprod(stp, grd))
    if (w$trace > 1) cat("Gradient projection = ",gprj)
##    tmp <- readline("   continue?")
    ## Do line search
    gvl<-w$lnsrch(fn,fbest, xb, stp, grv=grd, w, ...)
    print(str(gvl))
    if (attr(gvl, "FAIL")) {
        if (w$trace > 0) cat("Line search fails\n")
        break
    }
    fval <- attr(gvl,"Fval")
##    cat("after lnsrch w$nf now ",w$nf,"\n")
    if (w$trace > 1) {cat(" step =", gvl,"  fval=",fval ," w$nf=",w$nf,"\n")}
    xn<-xb+gvl*stp
    if (niter >= w$maxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
##    if (w$solver == "marquardt"){
##       if (fval <= fbest) {
##          xb <- xn
##          fbest <- fval
##          lambda <- lambda * w$lamdec
##          getH <- TRUE # ensure we start over
##       } else {
##          getH <- FALSE # don't want new H, grd
##          lambda <- lambda * w*laminc
##       }
##    }
    xb <- xn
    fbest <- fval
    niter <- niter + 1
    tmp <- readline("end iteration")   
  } # end repeat
  out<-NULL
  out$xs<-xn
  out$fv<-fval
  out$grd<-grd
  out$Hess<-H
  out$niter<-niter
  out 
}
