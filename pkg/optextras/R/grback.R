grback <- function(par, userfn, fbase=NULL, eps=1.0e-7, ...) {
   # Backward difference gradient approximation
   if (is.null(fbase)) fbase <- userfn(par, ...)  # ensure we function value at par
   df <- rep(NA, length(par))
   teps <- eps * (abs(par) + eps)
   for (i in 1:length(par)) {
      dx <- par
      dx[i] <- dx[i] - teps[i]
      df[i] <- (fbase - userfn(dx, ...))/teps[i]
   }
   df
}

