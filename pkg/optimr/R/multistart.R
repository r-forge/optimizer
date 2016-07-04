multistart <- function(parmat, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {
##
   nset <- nrow(parmat)
   npar <- ncol(parmat)
   if (nset < 1) stop("multistart: no starting parameters!")
   ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+7)
   ans.ret <- data.frame(ans.ret)
   pstring<-colnames(parmat)
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence", "kkt1", "kkt2", "xtimes")
   colnames(ans.ret)<-cnames
   row.names(ans.ret)<-1:nset

   for (imeth in 1:nset){
       start <- parmat[imeth, ]
       ans <- optimr(par=start, fn=fn, gr=gr, lower=lower, upper=upper, 
            method=method, hessian=hessian, control=control, ...)
       addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, 
                     ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
       
       ans.ret[imeth,] <- addvec
   }
   ans.ret

} ## end of multistart
