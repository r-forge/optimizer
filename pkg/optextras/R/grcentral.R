grcentral <- function(par, userfn, ...) {
   # Central difference gradient approximation
   eps <- .Machine$deriv.eps # DANGER: R core may use name
   df <- rep(NA, length(par))
   teps <- 0.5 * eps * (abs(par) + eps)
   for (i in 1:length(par)) {
      dp <- par
      dp[i]<-dp[i]+teps[i]
      dm <- par
      dm[i]<-dm[i]-teps[i]
      df[i] <- 0.5*(userfn(dp, ...) - userfn(dm,...))/teps[i]
   }
   df
}
 
