vmrf <- function(w,msetup=FALSE,...) {
## Fletcher 1970 variable metric style
   if (msetup) {
     w$c <- rep(0, w$npar) # old grad
     w$stp <- -1 # set step negative
#     w$H<-w$hess(w$xb,...)
#     w$nh <- w$nh + 1
     w$tdir <- w$c # on setup (??may not be needed)
     w$resetB <- TRUE
     w$grd<-w$gr(w$xb,...)
     #  if(w$trace > 1) {
     cat("vmrf setup - w$grd:")
     print(w$grd)
     # }
     w$ng <- w$ng + 1
     w$status <- "start"
     return(0)
   } else { # non setup ??still need to set gradient
     # BFGS update or reset
     if (w$resetB) {
       w$Bmat <- diag(w$npar) # set approx to inverse Hessian
       if (w$trace > 1) cat("Bmat reset to identity, no new gradient\n")  
       w$lastsd <- w$ng
     } else { # BFGS update of inverse Hessian approx Bmat
       w$c <- w$grd # "old" gradient
       w$grd<-w$gr(w$xb,...) # compute gradient
       cat("vmrf - w$grd:")
       print(w$grd)
       w$ng <- w$ng + 1
       cat("Check step and dirn exist: stp, tdir=",w$stp,"\n")
       print(w$tdir)
       w$tdir <- w$stp * w$tdir # ??add as.vector
       w$c <- w$grd - w$c # y
       D1 <- as.numeric(crossprod(w$tdir, w$c))
       if (w$trace > 1) cat("D1 =",D1," ")
       if (D1 <= 0) {
         if (w$trace > 1) cat("Update failed\n")
         w$Bmat <- diag(w$npar) # set approx to inverse Hessian
         if (w$trace > 1) cat("Bmat reset to identity\n")  
         w$lastsd <- w$ng
       } else {
         w$xx <- crossprod(w$Bmat, w$c)
         D2 <- as.numeric(crossprod(w$xx, w$c))
         D2 <- 1 + D2/D1 
#         cat("At update, D2=",D2," w$tdir:")
#         print((w$tdir))
#         cat("tdir * xx':\n")
#         print(w$tdir %*% t(w$xx))
#         cat("xx * tdir':\n")
#         print(w$xx %*% t(w$tdir))
#         cat("tdir * tdir':\n")
#         print(w$tdir %*% t(w$tdir))
         w$Bmat <- w$Bmat -
             (w$tdir %*% t(w$xx) + w$xx %*% t(w$tdir)- D2 * (w$tdir %*% t(w$tdir)))/D1
       } # end of matrix update
     } # end non-setup
     cat("At end update, Bmat:\n")
     print(w$Bmat)
     w$tdir <- - as.vector(crossprod(w$Bmat, w$grd))
     cat("tdir  vector? ",is.vector(w$tdir),":")
     print(w$tdir)
     tmp <- readline("?-")
     w$resetB <- FALSE # either way set so it doesn't do it again
   } # end of update
  ans <- w$tdir
  ans
}
