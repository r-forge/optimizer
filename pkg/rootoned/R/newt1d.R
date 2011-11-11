newt1d <- function(fn, gr, x0,
             offset = 1000.0, maxiter = 100, trace=FALSE,...){
#
# 1 dimensional Newton method for finding root xr of fn(x)
#  Need fn and gr(adient)
#
   epsilon <- .Machine$double.eps
   itn<-0
   xnew<-x0
   xold <- xnew+(1+abs(xnew)) # make sure it is different
   while ((xold + offset) != (xnew + offset)) {
       xold<-xnew
       fval<-fn(xold)
       gval<-gr(xold)
       itn<-itn+1
       if (itn > maxiter) stop("newt1d: Too many iterations")
       if (trace) cat(itn,":xold=",xold," f=",fval," g=",gval," xnew=",xnew,"\n")
       if (abs(gval) > (offset*epsilon*(abs(fval)+1.))) {
          xnew<-xold - fval/gval
       } else {
          stop("gradient too small")
       }
       if (trace) cat(itn,":xold=",xold," f=",fval," g=",gval," xnew=",xnew,"\n")
       ## cat("Change =",xnew-xold,"\n")
       ## tmp<-readline('next')
   }
   res<-list(root=xnew, froot=fval, itn=itn)    
}

