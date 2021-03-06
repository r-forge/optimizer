SNewton4<-function(x0,fn,gr,hess,control=c(trace=FALSE, maxit=1000),...) {
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

lnsrch<-function(fn,xc,d,...) { ## ?? replace this 
## Uses Brent's method to find the best stepsize gamma \in [0.1,1]
  flsch<-function(gm,fn,xc,d,...) {
    fval<-fn(xc+gm*d,...)
    fval
  }
  lout<-optimize(flsch,interval=c(0.1,1),fn=fn,xc=xc,d=d,...)
  rlout <- lout$min
  attr(rlout, "Fval")<- lout$objective
  rlout
}

  nmaxit<-2000 ## ?? replace
  eps0<-.Machine$double.eps
  eps <- 10*eps0
  lambda<-(eps0)^(1/4) ## ?? do better
  for (i in 1:nmaxit) { ## ?? replace with a while loop
    if (control$trace) {cat("Iteration ",i,":")}
    grd<-gr(x0,...)
    H<-hess(x0,...)
#    Haug<-H + (diag(H)+1.0)*lambda # To avoid singularity
#    stp<-solve(H, -grd) # In SNewton4 do this with svd
    ??
    ## Do line search
    gvl<-lnsrch(fn,x0,stp,...)
    if (control$trace) {cat(" step =", gvl,"  fval=", attr(gvl,"Fval"),"\n")}
    xn<-c(x0+gvl*stp)
    if ( max(abs(grd)) < eps ) break    
    if (i==nmaxit) {
      print("NewtonR: Failed to converge!")
      return(0)
    }
    x0<-xn
  }
  out<-NULL
  out$xs<-xn
  out$fv<-fn(xn,...)
  out$grd<-grd
  out$Hess<-H
  out$niter<-i
  out 
}




