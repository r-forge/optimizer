polyopt <- function(par, fn, gr=NULL, lower=NULL, upper=NULL, 
            methcontrol=NULL, hessian=FALSE, control=list(), ...) {
   print(str(methcontrol))
   nmeth <- nrow(methcontrol)
   cat(nmeth," Methods\n")
   npar <- length(par)
   if (nmeth < 1) stop("polyopt: no starting parameters!")
   ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+7)
   ans.ret <- data.frame(ans.ret)
   pstring<-names(par)
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence", "kkt1", "kkt2", "xtimes")
   colnames(ans.ret)<-cnames
   row.names(ans.ret)<-1:nmeth
   cpar <- par # use initial parameters here
# ?? may want a lot of checks here
   for (imeth in 1:nmeth){
       method <- methcontrol[[imeth,1]] # name of method
       cat("Method ",imeth," :",method,"\n")
       cat("imeth =",imeth,"\n")
       control$maxit <- methcontrol[[imeth,2]]
       control$maxfeval <- methcontrol[imeth,3]
       if (is.null(lower) && is.null(upper)) {
           cat("calling optimr unconstrained\n")
           ans <- optimr(cpar, fn=fn, gr=gr, 
            method=method, hessian=FALSE, control=control, ...)
       } else {
           cat("calling optimr bounded\n")
          ans <- optimr(cpar, fn=fn, gr=gr, lower=lower, upper=upper, 
            method=method, hessian=FALSE, control=control, ...)
       }
       ans$xtimes <- NA # 160703 -- not yet available
       addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, 
                     ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
       cpar <- ans$par # copy the parameters for next method
       ans.ret[imeth,] <- addvec
   }
   ans.ret

} ## end of polyopt
