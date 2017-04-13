snewton<-function(x0, fn, gr, hess, lower = NULL, upper = NULL, 
      control=list(lsmeth="default", solver="default", trace=2, maxit=500),...) {
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
  stepmax = 1.5,
  stepmin = -0.5,
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

halt <- function(curgval) {# test for termination
    halt <- FALSE # default is keep going
    # tests on too many counts
    if (niter > control$maxit) {
       if (control$trace > 0) cat("Too many (",niter," iterations\n")
       halt <- TRUE
       attr(halt,"convcode")<- 1
       return(halt)
    }
    
    if (nf > control$maxfevals){
       if (control$trace > 0) cat("Too many (",nf," function evaluations\n")
       halt <- TRUE
       attr(halt,"convcode")<- 91 # ?? value
       return(halt)
    }
#    if (ng > control$maxgevals){} # not implemented
#    if (nh > control$maxhevals){} # not implemented
    gmax <- max(abs(curgval))
    if (gmax <= control$epstol) {
      if (control$trace > 0) cat("Small gradient norm",gg,"\n")
      halt <- TRUE
      attr(halt,"convcode")<- 0 # OK
      return(halt)
    }
}


if (control$lsmeth == "default") {
  lnsrch<-function(fn, fbest, xc, d, grv, ...) { # Line search using internal optimize()
    cat("fn:\n")
    print(fn)
    ## Uses Brent's method to find the best stepsize in interval
    # fbest is best function value so far. NOT used.
    # grv is numeric gradient vector -- NOT used
    # ?? more documentation
    flsch<-function(st) {
      # computes the function value at stepsize st on line (xc + gm*d)
      # Essentially flsch(st)
      # gm: step size
      # fn: objective function
      # xc: base set of parameters
      # d : search direction
      nf <- nf +1
      fval<-fn(xc+st*d,...)
      fval
    }
    cat("function at ends of interval\n")
    sta <- control$stepmin
    cat("f(",sta,")=", flsch(sta),"\n")
    stb <- control$stepmax
    cat("f(",stb,")=", flsch(stb),"\n")
    
    lout<-optimize(flsch,interval=c(control$stepmin, control$stepmax),...)
    # ?? Need to count functions
    rlout <- lout$min
    attr(rlout, "Fval")<- lout$objective
    rlout # Note: returns stepsize, not x
  } # end default line search
} else if (control$lsmeth == "backtrack") {
  lnsrch<-function(fn, fbest, xc, d, grv, ...) { # backtrack line search
    # ?? count fevals?
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
   } # end backtrack line search
} else if (control$lsmeth == "none") { 
   lnsrch <- function(fn, fbest, xc, d, grv, ...) {
      rlout <- 1 # Does nothing! 
      attr(rlout, "Fval") <- fbest
      rlout
   }
}
  lambda<-control$lamstart ## ?? do better
  niter <- 1
  xb <- x0 # best so far
  fbest <- fn(xb, ...)
  newH <- TRUE
#  while (niter < control$maxit) { # main loop
  repeat {
    if (newH) {
        if (control$trace > 0) {cat("Iteration ",niter,":")}
        grd<-gr(xb,...)
        ng <- ng + 1
        H<-hess(xb,...)
        nh <- nh + 1
    }
#    if ( max(abs(grd)) < control$epstol ) break    
    if ( halt(grd) ) break
    if (control$solver == "default") {
      stp<-try(solve(H, -grd))
      if (class(stp) == "class-error") {
          stop("Failure of default solve of Newton equations")
      }
    } else if (control$solver == "marquardt") {
       Haug<-H + (diag(H)+1.0)*lambda # To avoid singularity
       stp <- solve(Haug, -grd)
    }
    if (control$trace > 0) {
         cat("Search vector:")
         print(stp)
    }
    gprj <- as.numeric(crossprod(stp, grd))
    cat("Gradient projection = ",gprj)
    tmp <- readline("   continue?")
    ## Do line search
    gvl<-lnsrch(fn,fbest, x0, stp, grv=NULL, ...)
# lnsrch<-function(fn, fbest, xc,d,grv, ...)
    fval <- attr(gvl,"Fval")
    if (control$trace > 0) {cat(" step =", gvl,"  fval=",fval ,"\n")}
    xn<-x0+gvl*stp
    if (niter >= control$maxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
    if (control$solver == "marquardt"){
       if (fval <= fbest) {
          xb <- xn
          fbest <- fval
          lambda <- lambda * control$lamdec
          newH <- TRUE # ensure we start over
       } else {
          newH <- FALSE # don't want new H, grd
          lambda <- lambda * control*laminc
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




