snewton<-function(x0, fn, gr, hess, control=list(trace=1, maxit=500),...) {
## Safeguarded Newton minimizer 
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

npar <- length(x0)
nf <- ng <- nh <- niter <- 0 # counters

ctrldefault <- list(
  trace = 0,
  maxit = 500,
  maxfevals = npar*500,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  stepdec = 0.2, 
  stepmax = 5,
  stepmin = 0,
  offset = 100.0,
  defstep=1,
  bigval = .Machine$double.xmax*0.01
)  

ncontrol <- names(control)
nctrld <- names(ctrldefault)
for (onename in nctrld) {
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrldefault[onename]
  }
}
trace <- control$trace # convenience

  xb <- x0 # best so far
  fbest <- fn(xb, ...)
  nf <- nf + 1 
  grd<-gr(xb,...) # compute gradient
  ng <- ng + 1
  fval <- control$bigval # to ensure comparison unfavourable
  #  while (niter < control$maxit) { # main loop
  repeat { # MAIN LOOP
    niter <- niter + 1
    if (trace > 1) cat("Termination test:")    
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter > control$maxit) {
      if (trace > 0) cat("Too many (",niter," iterations\n")
      halt <- TRUE
      convcode <- 1
      break
    }
    cat("tt nf=",nf,"\n")
    if (nf > control$maxfevals){
      if (trace > 0) cat("Too many (",nf," function evaluations\n")
      halt <- TRUE
      convcode <- 91 # ?? value
      break
    }
    #    if (ng > control$maxgevals){} # not implemented
    #    if (nh > control$maxhevals){} # not implemented
    gmax <- max(abs(grd))
    if (trace > 1) cat("current gradient norm =",gmax,"\n")
    if (gmax <= control$epstol) {
      if (trace > 1) cat("Small gradient norm\n")
      halt <- TRUE
      convcode <- 0 # OK
      break
    }
    # Note if we get here, 
    if (trace > 0) {cat("Iteration ",niter,":")}
    H<-hess(xb,...)
    nh <- nh + 1
    d<-try(solve(H, -grd))
    if (class(stp) == "class-error") {
          stop("Failure of default solve of Newton equations")
    }
    if (trace > 1) {
         cat("Search vector:")
         print(d)
    }
    gprj <- as.numeric(crossprod(d, grd))
    if (trace > 1) cat("Gradient projection = ",gprj)
    st <- control$defstep
    xnew <- xb + st*d # new point
    if ((control$offset+xnew) == (control$offset+xb)) {
        convcode <- 92 # no progress
        if (trace > 0) cat("No progress before linesearch!\n")
    }
    fval <- fn(xnew, ...)    
    nf <- nf + 1
    if (trace > 1) cat("f(xnew)=",fval,"\n")
    while ((fval > fbest + control$acctol*st*gprj) 
           && ((control$offset+xnew) != (control$offset+xb))) { # continue until satisfied
        st <- st * control$stepdec
        fval <- fn(xnew, ...)    
        nf <- nf + 1
        if (trace > 1) cat("* f(xnew)=",fval,"\n")
    }
    if ((control$offset+xnew) == (control$offset+xb)) {
        convcode <- 93 # no progress in linesearch
        if (trace > 0) cat("No progress IN linesearch!\n")
        break
    }
    if (trace > 1) cat("end major loop\n")  
    xb <- xnew
    fbest <- fval
    tmp <- readline("end iteration")   
  } # end repeat
  out<-NULL
  out$par<-xn
  out$value<-fval
  out$grad<-grd
  out$Hess<-H
  out$niter<-niter
  out$convcode <- convcode
  out 
}

