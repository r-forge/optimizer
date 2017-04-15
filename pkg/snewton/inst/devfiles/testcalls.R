# Try testing calls to see what is transferred (eventually test also ...)
# setup
 x0<-c(1,2,3,4)
 fnt <- function(x, fscale=10){
    yy <- length(x):1
    val <- sum((yy*x)^2)*fscale
 }
grt <- function(x, fscale=10){
    nn <- length(x)
    yy <- nn:1
#    gg <- rep(NA,nn)
    gg <- 2*(yy^2)*x*fscale
    gg
}

hesst <- function(x, fscale=10){
   nn <- length(x)
   yy <- nn:1
   hh <- diag(2*yy^2*fscale)
   hh
}

# Here want to set up ways to access other functions that are called inside testcalls



testcalls<-function(x0, fn, gr, hess, lower = NULL, upper = NULL, 
      control=list(lsmeth="default", solver="default", trace=2, maxit=500),...) {
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

  # This is purely a test routine to see what happens with different ways
  # to pass elements
  
npar <- length(x0)
nf <- ng <- nh <- niter <- 0 # counters

ctrldefault <- list(
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
  bigval = .Machine$double.xmax*0.01
)  

ncontrol <- names(control)
nctrld <- names(ctrldefault)
for (onename in nctrld) {
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrldefault[onename]
  }
}


  out 
}

