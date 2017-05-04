vmrf <- function(w,msetup=FALSE,...) {
## Fletcher 1970 variable metric style
   if (msetup) {
     w$grd<-w$gr(w$xb,...)
     cat("vmrf - w$grd:")
     print(w$grd)
     w$ng <- w$ng + 1
     w$c <- rep(0, w$npar) # old grad
     w$stp <- -1 # set step negative
#     w$H<-w$hess(w$xb,...)
#     w$nh <- w$nh + 1
     w$lastng <- 0 # last time we reset B -- force to identity
     w$tdir <- 0 # on setup
   } else { # non setup
     # BFGS update or reset
     if (w$lastng == 0) {
       w$Bmat <- diag(w$npar) # set approx to inverse Hessian
       if (w$trace > 1) cat("Bmat reset to identity\n")  
       w$lastng <- w$ng
     } else { # BFGS update of inverse Hessian approx Bmat
       D1 <- 0 
       w$tdir <- w$stp * w$tdir
       w$c <- w$grd - w$c # y
       D1 <- crossprod(w$tdir, w$c)
       if (D1 <= 0) {
         if (w$trace > 1) cat("Update failed\n")
         w$Bmat <- diag(w$npar) # set approx to inverse Hessian
         if (w$trace > 1) cat("Bmat reset to identity\n")  
         w$lastng <- w$ng
       } else {
         w$xx <- crossprod(w$Bmat, w$c)
         D2 <- as.numeric(crossprod(w$xx, w$c))
         D2 <- 1 + D2/D1         
         w$Bmat <- w$Bmat -
             (w$tdir %*% t(w$xx) + w$xx %*% t(w$tdir)- D2 * w$tdir %*% t(w$tdir))/D1
       } # end of matrix update
     } # end of update
     cat("At end update, Bmat:\n")
     print(Bmat)
     w$tdir <- - crossprod(w$Bmat, w$grd)
   } # end non-setup
  ans <- w$tdir
  ans
}
