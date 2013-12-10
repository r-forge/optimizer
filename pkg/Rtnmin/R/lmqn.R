lmqn <- function (x, sfun, maxit, maxfun, stepmx, accrcy, trace, ...) {
## ---------------------------------------------------------
##  truncated-newton method for unconstrained minimization
##  (customized version)
## ---------------------------------------------------------
#  global vectors hyk sk yk sr yr & scalars yksk yrsr
## ---------------------------------------------------------
##  set up
## ---------------------------------------------------------
## format compact
## format short e
  n<-length(x)
  # JN: Define globals here. Is it necessary to set up.
   gtn<-list(yrsr=0, yksk=0, yr = rep(0, n), yk = rep(0, n), sr = rep(0, n),  sk = rep(0, n))
   envjn<<-list2env(gtn)
# end globals
   eps <- .Machine$double.eps
   upd1 <- 1 
   ncg  <- 0 
   xnorm  <- max(abs(x)) # norm(x,'inf') 
   ierror <- 0 
   if (stepmx < sqrt(accrcy) || maxfun < 1) {
     ierror <- -1 
     xstar <- x  
     almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
              nfngr=ncg)
     return(almqn)
   }
## ---------------------------------------------------------
##  compute initial function value and related information
## ---------------------------------------------------------
   fg <- sfun(x, ...)
#%    print(fg)
   g<-fg$g
   f<-fg$f
   gnorm  <- max(abs(g)) ##  norm(g,'inf') 
   nf     <- 1 
   nit    <- 0 
   if (trace)  cat("Itn ",nit," ",nf," ",ncg, " ",f, " ", gnorm,"\n")
## ---------------------------------------------------------
##  check for small gradient at the starting point.
## ---------------------------------------------------------
   ftest <- 1 + abs(f) 
   if (gnorm < .01*sqrt(eps)*ftest) {
     ierror <- 0 
     xstar <- x 
     almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
     return(almqn)
   }
## ---------------------------------------------------------
##  set initial values to other parameters
## ---------------------------------------------------------
   n      <- length(x) 
   icycle <- n-1 
   toleps <- sqrt(accrcy) + sqrt(eps) 
   rtleps <- accrcy + eps 
   difnew <- 0 
   epsred <- .05 
   fkeep  <- f 
   conv   <- FALSE 
   ireset <- 0 
   ipivot <- 0 
## ---------------------------------------------------------
##  initialize diagonal preconditioner to the identity
## ---------------------------------------------------------
   d <- rep(1,n) # as a vector 
## ---------------------------------------------------------
##  ..........main iterative loop..........
## ---------------------------------------------------------
##  compute search direction
## ---------------------------------------------------------
   argvec <- c(accrcy, gnorm, xnorm) 
   mres  <- modlnp (d, x, g, maxit, upd1, ireset, bounds=FALSE, ipivot, argvec, sfun, ...) 
   p <- mres$p
   cat("p from first call to modlnp\n")
   print(p)
   tmp<-readline("cont.")

   gtp <- mres$gtp
   ncg1<-mres$ncg1
   d <- mres$dnew

#%    cat("d from modlnp:")
#%    print(d)
   ncg <- ncg + ncg1 
   while ( ! conv) { 
#%       cat("top while ")
#%       tmp <- readline("Top of iteration")
      oldg   <- g 
      pnorm  <- max(abs(p)) # norm(p,'inf') 
      oldf   <- f 
## ---------------------------------------------------------
##  line search
## ---------------------------------------------------------
      pe     <- pnorm + eps 
      spe    <- stepmx/pe 
#%       cat("gtp, spe:", gtp, spe,"\n")
      alpha0 <- step1 (f, gtp, spe) 
      reslin <- lin1 (p, x, f, alpha0, g, sfun, ...) 
##      [x, f, g, nf1, ierror, alpha] <- 
      x <- reslin$xnew # need fixup
      f <- reslin$fnew
      g <- reslin$gnew
      nf1 <- reslin$nf1
      ierror <- reslin$ierror
      alpha <- reslin$alpha1
      cat("after lin1, alpha=",alpha,"\n")
   tmp<-readline("cont.")

     nf <- nf + nf1 
## ---------------------------------------------------------
      nit <- nit + 1 
      gnorm <- max(abs(g)) # norm(g,'inf') 
### Display info
      if (trace) cat("Itn ",nit," ",nf," ",ncg, " ",f, " ", gnorm,"\n")
      if (ierror == 3) { 
         if (length(ncg) == 0) { ncg <- 0 } # ?? is.null(ncg)??
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
         return(almqn)
      }
## ---------------------------------------------------------
##  stop if more than maxfun evalutations have been made
## ---------------------------------------------------------
      if (nf >= maxfun) {
        ierror <- 2 
        xstar <- x 
        almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
        return(almqn)
      }
## ---------------------------------------------------------
##  set up for convergence and resetting tests
## ---------------------------------------------------------
      ftest  <- 1 + abs(f) 
      xnorm  <- max(abs(x)) # norm(x,'inf') 
      difold <- difnew 
      difnew <- oldf - f 
      envjn$yk     <- g - oldg 
      envjn$sk     <- alpha*p 
      if (icycle == 1) {
          if (difnew >   2*difold) { epsred <-   2*epsred }
          if (difnew < 0.5*difold) { epsred <- 0.5*epsred }
      }
## ---------------------------------------------------------
##  convergence test
## ---------------------------------------------------------
      conv <- ( ( (alpha*pnorm < toleps*(1 + xnorm)) &&
                 (abs(difnew) < rtleps*ftest)  &&
                 (gnorm < accrcy^(1/3)*ftest) )||
                 ( gnorm < .01*sqrt(accrcy)*ftest ) )
      if (conv) {
        ierror <- 0 
        xstar <- x 
        almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
        return(almqn)
      }
## ---------------------------------------------------------
##  update lmqn preconditioner
## ---------------------------------------------------------
      envjn$yksk <- as.numeric(crossprod(envjn$yk, envjn$sk) )
      ireset <- ((icycle == n-1) || (difnew < epsred*(fkeep-f)) )
      if ( ! ireset) { 
          envjn$yrsr <- as.numeric(crossprod(envjn$yr, envjn$sr))
          ireset <- (envjn$yrsr <= 0) 
      }
      upd1 <- (envjn$yksk <= 0) 
## ---------------------------------------------------------
##  compute search direction
## ---------------------------------------------------------
      argvec <- c(accrcy, gnorm, xnorm) 
##      cat("ireset =", ireset," upd1=",upd1,"\n")
#%       cat("New d to modlnp:")
#%       print(d)
##   [p, gtp, ncg1, d] <- ...
      mres <- modlnp (d, x, g, maxit, upd1, ireset, 0, ipivot, argvec, sfun, ...) 
      p <- mres$p
      gtp <- mres$gtp
      ncg1<-mres$ncg1
      d <- mres$dnew
      
#%       cat("returned dnew, envjn$john: ", envjn$john)
#%       print(d)
      ncg <- ncg + ncg1 
## ---------------------------------------------------------
##  store information for lmqn preconditioner
## ---------------------------------------------------------
      if (ireset) {
          envjn$sr <- envjn$sk 
          envjn$yr <- envjn$yk 
          fkeep <- f 
          icycle <- 1 
      } else {
          envjn$sr <- envjn$sr + envjn$sk 
          envjn$yr <- envjn$yr + envjn$yk 
          icycle <- icycle + 1 
      }
   } # end while
## [xstar, f, g, ierror] = ..
   almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
} # end lmqn

