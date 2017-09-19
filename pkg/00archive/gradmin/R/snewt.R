snewt <- function(w,msetup=FALSE,...) {
## Safeguarded Newton code for gradminu
   if (msetup) {
     w$grd<-w$gr(w$xb,...)
     if (w$trace > 2) {
        cat("snewt - w$grd:")
        print(w$grd)
     }
     w$ng <- w$ng + 1
     w$tdir <- 0 # on setup
   } else {
     # solve Newton equations
     w$H<-w$hess(w$xb,...)
     w$nh <- w$nh + 1
     w$tdir <- w$solver(w$H, -w$grd)
     if (class(w$tdir) == "class-error") {
       stop("Failure of default solve of Newton equations")
     }     
     w$grd<-w$gr(w$xb,...)
     if (w$trace > 2) {
       cat("snewt - w$grd:")
       print(w$grd)
     }
     w$ng <- w$ng + 1
   }
  ans <- w$tdir
  ans
}
