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
  
cat("control:")
print(control)
  
npar <- length(x0)

# set up workspace
ws <- list(
  lsmeth = "default",
  solver = "default",
  trace = 0,
  maxit = 500,
  maxfevals = npar*500,
  laminc = 10,
  lamdec = 0.4,
  lamstart = 0.01,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  svmin = 0.0,
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
  niter = 0
)  

ncontrol <- names(control)
nws <- names(ws)
for (onename in ncontrol) {
  if (onename %in% nws) {
    ws[onename]<-control[onename]
  }
}

cat("ws$trace=",ws$trace,"\n")

#cat("Workspace ws:")
#print(ws)
lnsrch <- lsback
## lnsrch <- lsnone # default to unit step
## if (ws$lsmeth == "lsback") {lnsrch <- lsback}
## else {stop("Undefined lsmeth = ",ws$lsmeth)}

## lnsrch <- pracma:fminbnd
## print(lnsrch)
## tmp <- readline("done")

## if (ws$lsmeth == "default") {
##   st <- lsback(fn, fbest, xc, d, grv, ...)
   ## } else if (ws$lsmeth == "backtrack") {
   ## } else if (ws$lsmeth == "none") { 
   ##    lnsrch <- function(fn, fbest, xc, d, grv, ...) {
   ##       rlout <- 1 # Does nothing! 
   ##       attr(rlout, "Fval") <- fbest
   ##       rlout
   ## }
   ## }
  lambda<-ws$lamstart ## ?? do better
  niter <- 1
  xb <- x0 # best so far
  fbest <- fn(xb, ...)
  ws$nf <- ws$nf + 1
  cat("ws$nf now ",ws$nf,"\n")
  newH <- TRUE
#  while (niter < ws$maxit) { # main loop
  repeat {
    if (newH) {
        if (ws$trace > 0) {cat("Iteration ",niter,":")}
        grd<-gr(xb,...)
        ws$ng <- ws$ng + 1
        H<-hess(xb,...)
        ws$nh <- ws$nh + 1
    }
    cat("Termination test:")    
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter >= ws$maxit) {
        if (ws$trace > 0) cat("Too many (",niter," iterations\n")
        halt <- TRUE
        convcode <- 1
        break
    }
    cat(" ws$nf=",ws$nf,"  ws$maxfevals=",ws$maxfevals,"\n")
    if (ws$nf >= ws$maxfevals) {
      cat("Stopping\n")
      if (ws$trace > 0) cat("Too many ",ws$nf," function evaluations\n")
      halt <- TRUE
      convcode <- 91 # ?? value
      break
    }
    #    if (ws$ng > ws$maxgevals){} # not implemented
    #    if (ws$nh > ws$maxhevals){} # not implemented
    gmax <- max(abs(grd))
    if (gmax <= ws$epstol) {
      if (ws$trace > 0) cat("Small gradient norm",gmax,"\n")
      halt <- TRUE
      convcode <- 0 # OK
      break
    }
    if (ws$solver == "default") {
      stp<-try(solve(H, -grd))
      if (class(stp) == "class-error") {
          stop("Failure of default solve of Newton equations")
      }
    } else if (ws$solver == "marquardt") {
       Haug<-H + (diag(H)+1.0)*lambda # To avoid singularity
       stp <- solve(Haug, -grd)
    }
    if (ws$trace > 0) {
         cat("Search vector:")
         print(stp)
    }
    gprj <- as.numeric(crossprod(stp, grd))
    cat("Gradient projection = ",gprj)
    tmp <- readline("   continue?")
    ## Do line search
    gvl<-lnsrch(fn,fbest, xb, stp, grv=grd, ws, ...)
# lnsrch<-function(fn, fbest, xc,d,grv, ...)
    print(str(gvl))
    fval <- attr(gvl,"Fval")
    cat("after lnsrch ws$nf now ",ws$nf,"\n")
    if (ws$trace > 0) {cat(" step =", gvl,"  fval=",fval ," ws$nf=",ws$nf,"\n")}
    xn<-xb+gvl*stp
    if (niter >= ws$maxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
    if (ws$solver == "marquardt"){
       if (fval <= fbest) {
          xb <- xn
          fbest <- fval
          lambda <- lambda * ws$lamdec
          newH <- TRUE # ensure we start over
       } else {
          newH <- FALSE # don't want new H, grd
          lambda <- lambda * ws*laminc
       }
    }
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
