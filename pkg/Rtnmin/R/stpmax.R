stpmax <- function(stepmx, pe, x, p, ipivot, low, up) {
##------------------------------------------------
## compute the maximum allowable step length
## (spe is the standard (unconstrained) max step)
##------------------------------------------------
   ## cat("stpmax: stepmx, pe:",stepmx, pe,"\n")

   spe  <- stepmx / pe 
   al   <- spe 
   au   <- spe 
##------------------------------------------------
   indl <- which(ipivot==0 & p < 0) 
   if ( length(indl)>0) {
      tl   <- low[indl] - x[indl] 
      al   <- min(tl/p[indl]) 
   }
##------------------------------------------------
   indu <- which(ipivot==0 & p > 0) 
   if (length(indu)>0) {
      tu   <- up[indu] - x[indu] 
      au   <- min(tu/p[indu]) 
   }
##------------------------------------------------
   spe  <- min(c(spe, al, au)) 
}