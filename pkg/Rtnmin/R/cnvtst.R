cnvtst <- function  (alpha, pnorm, xnorm, 
		 dif, ftest, gnorm, gtp, f, flast, g, ipivot, accrcy) {

## ---------------------------------------------------------
##  test for convergence
## ---------------------------------------------------------
##  set up
## ---------------------------------------------------------
conv   <- 0;
eps <- .Machine$double.eps
toleps <- sqrt(accrcy) + sqrt(eps);
rtleps <- accrcy + eps;
imax   <- 0;
ltest  <- (flast - f <= -0.5*gtp);
## ---------------------------------------------------------
##  if anti-zigzag test satisfied, test multipliers;
##  if appropriate, modify active set
## ---------------------------------------------------------
   if ( ! ltest) {
      ind   <- which( (ipivot != 0)  & (ipivot != 2))
      if ( length(ind) > 0 ) {  ## how to ensure ind not empty
         t <- -sum(ipivot[ind]*g[ind])
         cmax <- min(t) # ?? why [cmax, imax] = min(t) ??
         imax <- cmax
         if (cmax >= 0) { imax <- 0 }
      }
   }
   if (imax != 0) {
      ipivot[ ind[imax] ] <- 0;
      flast <- f;
   } else {
      conv = ( ( (alpha*pnorm < toleps*(1 + xnorm) )
           && (abs(dif) < rtleps*ftest)
           && (gnorm < accrcy^(1/3)*ftest) ) || 
            ( gnorm < .01*sqrt(accrcy)*ftest) )
   }
   result <- list(conv=conv, flast1=flast, ipivot1=ipivot)
} 
