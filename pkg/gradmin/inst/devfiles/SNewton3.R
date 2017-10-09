SNewton3<-function(x0,fn,gr,hess,control=c(trace=FALSE, maxit=1000),...) {
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

  nmaxit<-2000 ## ?? replace
  eps0<-.Machine$double.eps
  eps <- 10*eps0
  n <- length(x0) # number of parameters
  lambda<-(eps0)^(1/4) ## ?? do better
  phi <- 1.0
  laminc <- 10
  lamdec <- 0.4
  f0 <- fn(x0, ...)
  for (itn in 1:nmaxit) { ## ?? replace with a while loop
    if (control$trace) {cat("Iteration ",itn,":")}
    grd<-gr(x0,...)
    if ( max(abs(grd)) < eps ) break # ?? change in SNewton2 as well
    H<-hess(x0,...)
    repeat { # loop until lower function found
#      Haug<-H + (diag(H)+phi*diag(n))*lambda # To avoid singularity
      Haug<-H + (diag(n))*lambda # To avoid singularity
      #    stp<-solve(H, -grd) # SNewton2 approach
# Here we will use Marquardt
      stp<-solve(Haug, -grd) # SNewton2 approach
      xn<-c(x0+stp)
      if (itn >= nmaxit) {
        print("NewtonR: Failed to converge!")
        return(0)
      }
      fval <- fn(xn, ...)
      cat("    lambda = ",lambda,"  fval =",fval,"\n")
      if (fval < f0) {
         x0<-xn
         f0 <- fval
         lambda <- lamdec * lambda
         break
      } else {
         lambda <- lambda * laminc
      }
    } # end repeat
  } 
  out<-NULL
  out$xs<-x0
  out$fv<-f0
  out$grd<-grd
  out$Hess<-H
  out$niter<-itn
  out 
}
