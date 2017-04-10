snewton<-function(x0, fn, gr, hess, lower = NULL, upper = NULL, 
      control=list(lsmeth="default", solver="default", trace=FALSE, maxit=500),...) {
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
nf <- ng <- nh <- 0 # counters

ctrldefault <- list(
  lsmeth = "default",
  solver = "default",
  trace = FALSE,
  maxit = 500,
  maxfevals = npar*500,
  laminc = 10,
  lamdec = 0.4,
  lamstart = 0.01,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  svmin = 0.0,
  stepdec = 0.2, 
  stepmax = 1.2,
  stepmin = 0.1,
  offset = 100.0,
  bigval = .Machine$double.xmax*0.01
)  

ncontrol <- names(control)
nctrld <- names(ctrldefault)
for (onename in nctrld) {
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrldefault[onename]
  }
}

flsch<-function(gm,fn,xc,d,...) {
  # computes the function value at stepsize gm on line (xc + gm*d)
  # gm: step size
  # fn: objective function
  # xc: base set of parameters
  # d : search direction
  fval<-fn(xc+gm*d,...)
  nf <- nf +1
  fval
}

if (control$lsmeth == "default") {
    lnsrch<-function(fn, fbest, xc,d,grv, ...) { # Line search using internal optimize()
      ## Uses Brent's method to find the best stepsize gamma \in [0.1,1]
      # fbest is best function value so far. NOT used.
      # grv is numeric gradient vector -- NOT used
      # ?? more documentation
    lout<-optimize(flsch,interval=c(stepmin, stepmax),fn=fn,xc=xc,d=d,...)
    rlout <- lout$min
    attr(rlout, "Fval")<- lout$objective
    rlout # Note: returns stepsize, not x
  } # end default line search
} else if (lsmeth == "backtrack") {
  lnsrch<-function(fn, fbest, xc, d, grv, ...) { # backtrack line search
    st <- 1.0
    gproj <- as.numeric(crossprod(grv, xc) )
    repeat {
      xnew <- xc + st*d # new point
      if ((offset+xnew) == (offset+xc)) { # no better parameters
          st <- 0
          rlout <- st
          attr(rlout,"Fval")<-fbest # Assume we pass this in
          return(rlout)
      }
      fval <- flsch(xnew, ...)
      if (control$trace > 1) cat("Step = ",st," fval = ",fval,"\n")
      if (fval <= fbest + acctol*st*gproj) break
      st <- stepdec*st # new step
    }
    rlout <- st
    attr(rlout, "Fval")<- fval
    rlout
   } # end default line search

}

  nmaxit<-2000 ## ?? replace
  eps0<-.Machine$double.eps
  eps <- 10*eps0
  lambda<-(eps0)^(1/4) ## ?? do better
  for (itn in 1:nmaxit) { ## ?? replace with a while loop
    if (control$trace) {cat("Iteration ",itn,":")}
    grd<-gr(x0,...)
    if ( max(abs(grd)) < eps ) break    
    H<-hess(x0,...)
#    Haug<-H + (diag(H)+1.0)*lambda # To avoid singularity
    stp<-solve(H, -grd)
    ## Do line search
    gvl<-lnsrch(fn,x0,stp,...)
    fval <- attr(gvl,"Fval")
    if (control$trace) {cat(" step =", gvl,"  fval=",fval ,"\n")}
    xn<-c(x0+gvl*stp)
    if (itn >= nmaxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
    x0 <- xn
    f0 <- fval
  }
  out<-NULL
  out$xs<-xn
  out$fv<-fval
  out$grd<-grd
  out$Hess<-H
  out$niter<-itn
  out 
}




