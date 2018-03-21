#Example 2: Wood function
#
wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  return(res)
}
#gradient:
wood.g <- function(x){
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}
#hessian:
wood.h <- function(x){
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}

wood.fgh <- function(x){
      fval <- wood.f(x)
      gval <- wood.g(x)
      hval <- wood.h(x)
      attr(fval,"gradient") <- gval
      attr(fval,"hessian")<- hval
      fval
}

newton<-function(par, fn, gr, hess, control=list(trace=1, maxit=500),...) {
  ## Safeguarded Newton minimizer 
  ##
  ##Input
  ##       - fn is the function we wish to minimize
  ##?? fixup documentation here??
  ##       - par is the initial value
  ##       - ... is data used in the function fn
  ##Output (list) -- need to match optim() output!! ???
  ##       - xs is the value at the minimum
  ##       - fv is the fn evaluated at xs
  ##       - grd is the gradient
  ##       - Hess is the Hessian
  ##       - niter is the number of interations needed (gradient and Hessian evals).
  ##       - add fevals??, other reports
  
  npar <- length(par)
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
  cat("trace =",trace,"\n")
  
  xb <- par # best so far
  fbest <- fn(xb, ...)
  nf <- nf + 1 
  if (trace > 0) cat("Initial function value = ",fbest,"\n")
  if (trace > 1) print(xb)
  #  fval <- control$bigval # to ensure comparison unfavourable
  #  while (niter < control$maxit) { # main loop
  repeat { # MAIN LOOP
    niter <- niter + 1
    grd<-gr(xb,...) # compute gradient
    ng <- ng + 1
    if (trace > 1) cat("Termination test:")    
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter > control$maxit) {
      if (trace > 0) cat("Too many (",niter," iterations\n")
      halt <- TRUE
      convcode <- 1
      break
    }
    cat("nf=",nf,"\n")
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
    if (class(d) == "class-error") {
      stop("Failure of default solve of Newton equations")
    }
    if (trace > 1) {
      cat("Search vector:")
      print(d)
    }
    gprj <- as.numeric(crossprod(d, grd))
    if (trace > 1) cat("Gradient projection = ",gprj,"\n")
    st <- control$defstep
    xnew <- xb + st*d # new point
    if (all((control$offset+xnew) == (control$offset+xb))) {
      convcode <- 92 # no progress
      if (trace > 0) cat("No progress before linesearch!\n")
    }
    fval <- fn(xnew, ...)    
    nf <- nf + 1
    if (trace > 1) cat("f(xnew)=",fval,"\n")
    if (trace > 1) cat("end major loop\n")  
    xb <- xnew
    fbest <- fval
    tmp <- readline("end iteration")   
  } # end repeat
  out<-NULL
  out$par<-xb
  out$value<-fbest
  out$grad<-grd
  out$Hess<-H
  out$niter<-niter
  out$convcode <- convcode
  out
}


 
#################################################
x0 <- c(-3,-1,-3,-1) # Wood standard start

library(snewton)
cat("This FAILS to find minimum\n")
wd <- newton(x0, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=2))
print(wd)

cat("\n\n nlm() gives same result\n")
tnlm <- nlm(wood.fgh, x0, print.level=2)
print(tnlm)


