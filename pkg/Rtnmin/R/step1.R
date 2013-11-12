step1 <- function(f, gtp, smax) {
##---------------------------------------------------------
## step1 returns the length of the initial step to be 
## taken along the vector p in the next linear search.
##---------------------------------------------------------
## [fm is supposed to be an estimate of the optimal function value]
   eps<-.Machine$double.eps
   fm <- 0
   d  <- abs(f-fm)
   alpha <- min(1, smax)
   if ((2*d <= (-gtp)) && (d >= eps)) {
      alpha = min(-2*d/gtp, smax)
   }
   alpha # need to ensure returned value
}