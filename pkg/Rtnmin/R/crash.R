crash <- function (x, low, up) {
##---------------------------------------------------------
## this initializes the constraint information, and
## ensures that the initial point satisfies 
##      low <= x <= up.
## the constraints are checked for consistency.
##---------------------------------------------------------
ierror <- 0 
if (any(low > up)) { ierror = - max(which(low > up))  } 
   # above is check on error in bounds specification
   xnew  <- pmax (low, x) # force params into bounds
   xnew  <- pmin (up, xnew) # No diagnostic! ??
   # we output revised parameters
   n <- length(x)
   ind <- which(low == x) 
   ipivot <- rep(0, n) ## zeros(size(x)) in Matlab
   if (length(ind) > 0) { ipivot[ind] <- -1 }
   ind  <- which(x == up) 
   if (length(ind) > 0) { ipivot[ind] <-  1 }
   ind <- which(low == up) 
   if (length(ind) > 0) { ipivot[ind] <- 2 }
   list(ipivot=ipivot, ierror=ierror, xnew=xnew)
}


