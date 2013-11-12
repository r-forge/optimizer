modz <- function (x, p, ipivot, low, up, flast, f, alpha){
##---------------------------------------------------------------------
## update the constraint matrix if a new constraint is encountered
##---------------------------------------------------------------------
   eps <- .Machine$double.eps
   indl <- which(ipivot == 0 & p < 0)
   if (length(indl) > 0) {
      toll <- 10 * eps * (abs(low[indl]) + 1)
      hitl <- which(x[indl]-low[indl] <= toll)
      if (length(hitl) > 0) {
         flast <- f
         ipivot[indl[hitl]] <- -1
      }
   }
##---------------------------------------------------------------------
   indu <- which((ipivot == 0) & (p > 0));
   if (length(indu) > 0) {
      tolu <- 10 * eps * (abs( up[indu]) + 1)
      hitu <- which(up[indu]-x[indu]  <= tolu)
      if (length(hitu) > 0) {
         flast <- f
         ipivot[indu[hitu]] <- 1
      }
   }
##---------------------------------------------------------------------
   flast1  <- flast
   ipivot1 <- ipivot
   list(flast1 = flast, ipivot1=ipivot)
}

