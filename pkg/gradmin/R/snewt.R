snewt <- function(w,msetup=FALSE,...) {
## Safeguarded Newton code for gradminu
   if (msetup) {
     w$grd<-w$gr(w$xb,...)
     cat("snewt - w$grd:")
     print(w$grd)
     w$ng <- w$ng + 1
     w$H<-w$hess(w$xb,...)
     w$nh <- w$nh + 1
     w$tdir <- 0 # on setup
   } else {
     # solve Newton equations
     w$tdir <- w$solver(w$H, -w$grd)
     if (class(w$tdir) == "class-error") {
       stop("Failure of default solve of Newton equations")
     }     
   }
  ans <- w$tdir
  ans
}
