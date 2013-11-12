gtims <- function (v, x, g, accrcy, xnorm, sfun, ...) {
##---------------------------------------------------------
## compute the product of the Hessian times the vector v;
## store result in the vector gv 
## (finite-difference version)
##---------------------------------------------------------
##   cat("v, x, g:")
##   print(v)
##   print(x)
##   print(g)
   delta <- sqrt(accrcy)*(1 + xnorm)/sqrt(sum(v^2))
   hg <- x + delta*v
   tresult<-sfun(hg, ...)
##   cat("tresult: ")
##   print(tresult)
#   gv <- sfun(hg, ...)$gval
   gv <- tresult$g
   gv <- (gv - g)/delta
##   cat("returning gv:")
##   print(gv)
   gv 
}
