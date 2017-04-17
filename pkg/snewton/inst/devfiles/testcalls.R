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

Newtdir <- function(x, fn, gr, hess, wspace, control=list(), ...) {
 # wspace <- control$wrk # This is more certain, but not much more effort
#  wrkenv <- parent.env(environment()) # the env of calling routine
  cat("Newtdir -- wspace:")
  print(wspace)
  cat("Newtdir -- fn, gr, hess:")
  print(fn)
  print(gr)
  print(hess)
  cat("Newtdir -- wspace$gcurr:")
  print(wspace$gcurr)
#  cat('exists("wspace$gcurr") =', exists("wspace$gcurr"),"\n")
  # seems exists is only for the NAME
  cat("is.null('wspace$gcurr') =", is.null('wspace$gcurr'),"\n")
  # Try to compute Newton direction
  if (is.null(wspace$gcurr)) {
     cat("Compute gradient gcurr in Newtdir\n")
     gcurr <- gr(x, ...)
     wspace$gcurr <- gcurr
  } else {
    cat("gcurr from wspace:")
    gcurr <- wspace$gcurr
    print(gcurr)
  }
  # Need some way to pass this back to the calling environment
  if (is.null(wspace$hcurr)) {
    hcurr <- hess(x, ...) # ditto
    wspace$hcurr <- hcurr
    cat("generated hcurr:")
    print(hcurr)
  } else { hcurr <- wspace$hcurr }
  cat("Newtdir, hcurr and gcurr\n")
  print(hcurr)
  print(gcurr)
  sdir <- solve(hcurr, -gcurr)
  attr(sdir,"wspace") <- wspace
  print(sdir)
  sdir
}
  


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

wspace <<- list2env(control) # NOTE: Global!

cat("generate gcurr in testcalls\n")
x <- x0
wspace$gcurr <- gr(x0, ...)
print(wspace$gcurr)

sd1 <- Newtdir(x, fn, gr, hess, wspace, control, ...)
cat("sd1:")
print(sd1)

cat("Now NULL gcurr in wspace\n")
wspace$gcurr <- NULL

sd2 <- Newtdir(x, fn, gr, hess,  wspace, control, ...)
cat("sd2:")
print(sd2)
cat("After sd2, wspace$gcurr:")
print(wspace$gcurr)

out <- sd2
out
}

cat("Testcalls -- 170415 -- to illustrate how R transfers information\n")

res10 <- testcalls(x0, fnt, grt, hesst, fscale=10.0)
cat("res10:   ")
print(res10)

cat("\nNow try with fscale=1\n")
res1 <- testcalls(x0, fnt, grt, hesst, fscale=1.0)
cat("res1:    ")
print(res1)

