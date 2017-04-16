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
  
  
npar <- length(x0)
nf <- ng <- nh <- niter <- 0 # counters

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
  bigval = .Machine$double.xmax*0.01
)  

ncontrol <- names(control)
nws <- names(ws)
for (onename in ncontrol) {
  if (! (onename %in% nws)) {
    ws[onename]<-control[onename]
  }
}

cat("Workspace ws:")
print(ws)

lnsrch <- pracma::fminbnd
print(lnsrch)
tmp <- readline("done")

if (ws$lsmeth == "default") {
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
#      nf <- nf +1
      fval<-fn(xc+st*d,...)
      fval
    }
    cat("function at ends of interval\n")
    sta <- ws$stepmin
    cat("f(",sta,")=", flsch(sta),"\n")
    stb <- ws$stepmax
    cat("f(",stb,")=", flsch(stb),"\n")
    
    #  lout<-optimize(flsch,interval=c(ws$stepmin, ws$stepmax),
    #                  lower=ws$stepmin, upper=ws$stepmax,...)
    # note fmin rather than objective in return  
    lout<-pracma::fminbnd(flsch,ws$stepmin, ws$stepmax, ...)
    cat("lnsrch lout:")
    print(lout)
    rlout <- lout$xmin
    #  cat("structure of rlout")
    #  print(str(rlout))
    attr(rlout, "Fval") <- lout$fmin
    attr(rlout, "fcount") <- (lout$niter + 1) # fevals is iterations + 1
    rlout # Note: returns stepsize, not x
  } # end default line search
} else if (ws$lsmeth == "backtrack") {
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
      if (ws$trace > 1) cat("Step = ",st," fval = ",fval,"\n")
      if (fval <= fbest + acctol*st*gproj) break
      st <- stepdec*st # new step
    }
    rlout <- st
    attr(rlout, "Fval")<- fval
    rlout
   } # end backtrack line search
} else if (ws$lsmeth == "none") { 
   lnsrch <- function(fn, fbest, xc, d, grv, ...) {
      rlout <- 1 # Does nothing! 
      attr(rlout, "Fval") <- fbest
      rlout
   }
}
  lambda<-ws$lamstart ## ?? do better
  niter <- 1
  xb <- x0 # best so far
  fbest <- fn(xb, ...)
  nf <- nf + 1
  newH <- TRUE
#  while (niter < ws$maxit) { # main loop
  repeat {
    if (newH) {
        if (ws$trace > 0) {cat("Iteration ",niter,":")}
        grd<-gr(xb,...)
        ng <- ng + 1
        H<-hess(xb,...)
        nh <- nh + 1
    }
    cat("Termination test:")    
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter > ws$maxit) {
        if (ws$trace > 0) cat("Too many (",niter," iterations\n")
        halt <- TRUE
        convcode <- 1
        break
    }
    cat("tt nf=",nf,"\n")
    if (nf > ws$maxfevals){
    if (ws$trace > 0) cat("Too many (",nf," function evaluations\n")
        halt <- TRUE
        convcode <- 91 # ?? value
        break
    }
    #    if (ng > ws$maxgevals){} # not implemented
    #    if (nh > ws$maxhevals){} # not implemented
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
    gvl<-lnsrch(fn,fbest, xb, stp, grv=NULL, ...)
# lnsrch<-function(fn, fbest, xc,d,grv, ...)
    print(str(gvl))
    fval <- attr(gvl,"Fval")
    nf <- nf + attr(gvl, "fcount")
    if (ws$trace > 0) {cat(" step =", gvl,"  fval=",fval ," nf=",nf,"\n")}
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

t1 <- gradminu(x0, fnt, grt, hesst, control=list(trace=2), fscale=3.0)
t1
