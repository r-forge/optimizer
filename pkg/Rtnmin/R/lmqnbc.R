lmqnbc <- function (x, sfun, lower, upper, maxit, maxfun, stepmx, accrcy, trace, ...) {
## ---------------------------------------------------------
##  This is a bounds-constrained truncated-newton method.
##  The truncated-newton method is preconditioned by a 
##  limited-memory quasi-newton method (computed by
##  this routine) with a diagonal scaling (routine ndia3).
##  For further details, see routine tnbc.
## ---------------------------------------------------------
##  global hyk sk yk sr yr yksk yrsr
## ---------------------------------------------------------
##  check that initial x is feasible and that the bounds 
##  are consistent
## ---------------------------------------------------------
   n <- length(x)

# JN: Define globals here
   gtn<-list(yrsr=0, yksk=0, yr = rep(0, n), yk = rep(0, n), 
        sr = rep(0, n),  sk = rep(0, n),
        hg=rep(0,n), hyk=rep(0,n), hyr=rep(0,n) )
   envjn<<-list2env(gtn)
# end globals

##   [ipivot, ierror, x] = crash(x, lower, upper) 
   cat("lmqnbc -- crout:")
   crout<-crash(x, lower, upper)
   ## print(crout)
   ierror <- crout$ierror
   ipivot <- crout$ipivot
   x<-crout$xnew # in case x changed by bounds
   f = 0 
   g = rep(0,n) 
   if (ierror != 0) {
      stop('LMQNBC: terminating (no feasible point)')
      ##   return 
   }
## ---------------------------------------------------------
##  initialize variables, parameters, and constants
## ---------------------------------------------------------
   options(digits=5)
   ## cat("stempmx, accrcy, maxfun:",stepmx, accrcy, maxfun, "\n")
   cat('  it     nf     cg           f             |g|\n')
   eps <- .Machine$double.eps
   upd1   <- TRUE
   ncg    <- 0 
   conv   <- FALSE
   xnorm  <- max(abs(x))
   ierror <- 0 

   if ( (stepmx < sqrt(accrcy)) || (maxfun < 1) ) { 
      ## cat("stepmx < sqrt(accrcy) terminator\n")
      ierror <- -1 
      xstar <- x  
      almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                 nfngr=ncg)
     cat("Exiting lmqnbc - almgn:")
     print(almqn)
     return(almqn)  # return 
   }
## ---------------------------------------------------------
##  compute initial function value and related information
## ---------------------------------------------------------
#   cat("Try initial fn\n")
   fg<- sfun(x, ...)
   nf     <- 1 
   nit    <- 0 
   g<-fg$g
   f<-fg$f
   flast  <- f
   gnorm  <- max(abs(g)) ##  norm(g,'inf') 
## ---------------------------------------------------------
##  Test if Lagrange multipliers are non-negative.
##  Because the constraints are only bounds, the Lagrange
##  multipliers are components of the gradient.
##  Then form the projected gradient.
## ---------------------------------------------------------
   ind <- which((ipivot != 2) & 
             (as.numeric(crossprod(ipivot,g)) >0 ) ) 
   if (length(ind) > 0) {  
      ipivot[ind] <- rep(0, length(ind)) 
   } 
   g <- ztime (g, ipivot) 
   gnorm <- max(abs(g)) 
   cat(nit,"\t", nf,"\t", ncg,"\t", f,"  ", gnorm,"\n")
   ## print(x)
   ## print(g)
## ---------------------------------------------------------
##  check if the initial point is a local minimum.
## ---------------------------------------------------------
   ftest <- 1 + abs(f) 
   if (gnorm < .01*sqrt(eps)*ftest) {
      xstar <- x 
      almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
              nfngr=ncg)
      return(almqn)
   }

      ## cat("before set initial, flast=",flast,"\n")

## ---------------------------------------------------------
##  set initial values to other parameters
## ---------------------------------------------------------
   icycle <- n-1 
   ireset <- 0 
   bounds <- TRUE 
   difnew <- 0 
   epsred <- .05 
   fkeep  <- f 
   d      <- rep(1,n) 
## ---------------------------------------------------------
##  ..........main iterative loop..........
## ---------------------------------------------------------
##  compute the search direction
## ---------------------------------------------------------
   argvec <- c(accrcy, gnorm, xnorm) 
##[p, gtp, ncg1, d] <- ...
##	modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, ##   argvec, sfun)
   mres  <- modlnp (d, x, g, maxit, upd1, ireset,
             bounds, ipivot, argvec, sfun, ...) 
   ncg1 <- mres$ncg1
   gtp  <- mres$gtp
   d    <- mres$dnew
   p    <- mres$p
   ## cat("p:")
   ## print(p)
   ## tmp<-readline("cont.")

   ncg <- ncg + ncg1 
   while (!conv) {
      oldg <- g 
      pnorm <- max(abs(p)) # norm2(p, 'inf') 
      oldf <- f 
## ---------------------------------------------------------
##  line search
## ---------------------------------------------------------
      pe <- pnorm + eps 
      spe <- stpmax (stepmx, pe, x, p, ipivot, lower, upper) 
      ## cat("spe=",spe,"\n")
      alpha <- step1 (f, gtp, spe) 
      alpha0 <- alpha 
      ## cat("alpha0=",alpha0,"\n")
##   [x_new, f_new, g_new, nf1, ierror, alpha] <- lin1 (p, x,
##       f, alpha0, g, sfun) 
      reslin <- lin1 (p, x, f, alpha0, g, sfun, ...) 
      ierror <- reslin$ierror
      alpha <- reslin$alpha1
## ---------------------------------------------------------
      if ((alpha == 0) && (alpha0 != 0) || (ierror == 3)){ 
          cat('Error in Line Search\n') 
          cat('    ierror = ', ierror, "\n") 
          cat('    alpha  = ',alpha, "\n") 
          cat('    alpha0 = ', alpha0, "\n") 
          cat('    gtp    = ', gtp, "\n") 
          ## ############################
          cat('    |g|     = ', norm2(g), "\n") 
          cat('    |p|     = ', norm2(p), "\n") 
          tmp <- readline('Hit any key to continue')
          ## ############################
      } 
      ## #######################
      x <- reslin$xnew # need fixup
      f <- reslin$fnew
      g <- reslin$gnew
      nf1 <- reslin$nf1
      nf  <- nf  + nf1 
      nit <- nit +   1 
## ---------------------------------------------------------
##  update active set, if appropriate
## ---------------------------------------------------------
      newcon <- FALSE 
      ## cat("Check active set - alpha, spe, f:",alpha, spe,f,"\n")
      if (abs(alpha-spe) <= 10*eps) {
         newcon <- TRUE
         ierror <- 0 
          ## cat("flast:",flast,"\n")
#         [ipivot, flast] <- modz (x, p, ipivot, lower, upper, 
#             flast, f, alpha) 
         ## cat("ipivot before and after update by modz:\n")
         ## print(ipivot)
         modzres<-modz(x, p, ipivot, lower, upper, 
                      flast, f, alpha) 
         ipivot <- modzres$ipivot1
         flast<-modzres$flast1
         ## print(ipivot)
      }
      if (ierror == 3) { 
         xstar <- x  
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                 nfngr=ncg)
         cat("Error updating active set -- almqn:")
         print(almqn)
         return(almqn)
      }
## ---------------------------------------------------------
##  stop if more than maxfun evaluations have been made
## ---------------------------------------------------------
      if (nf > maxfun) { 
         ierror <- 2 
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                     nfngr=ncg)
         cat("Too many function evaluations -- almqn:")
         print(almqn)
         return(almqn)
      } 
## ---------------------------------------------------------
##  set up for convergence and resetting tests
## ---------------------------------------------------------
      difold <- difnew 
      difnew <- oldf - f # scalar 
      if (icycle == 1) {
         if (difnew >  2*difold) { epsred <-  2*epsred } 
         if (difnew < .5*difold) { epsred <- .5*epsred } 
      } 
      gv    <- ztime (g, ipivot) 
      gnorm <- max(abs(gv)) 
      ftest <- 1 + abs(f) 
      xnorm <- max(abs(x))
############### DISPLAY ############## 
      cat(nit,"\t", nf,"\t", ncg,"\t", f,"  ", gnorm,"\n")
      ## print(x)
      ## print(g)
## ---------------------------------------------------------
##  test for convergence
## ---------------------------------------------------------
##   [conv, flast, ipivot] <- cnvtst (alpha, pnorm, xnorm, ...
##	    difnew, ftest, gnorm, gtp, f, flast, g, ...
##	    ipivot, accrcy) 
      ## cat("before cnvtst, flast=",flast,"\n")
      ctres <- cnvtst (alpha, pnorm, xnorm, 
	       difnew, ftest, gnorm, gtp, f, flast, g, 
	       ipivot, accrcy) 
      conv <- ctres$conv
      flast <- ctres$flast1
      ipivot <- ctres$ipivot1
      ## cat("after cnvtst - conv, flast, ipivot:", conv, flast,"\n")
      ## print(ipivot)
      if (conv) {
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                     nfngr=ncg)
         return(almqn)
      } 
      g <- ztime (g, ipivot) 
## ---------------------------------------------------------
##  modify data for LMQN preconditioner
## ---------------------------------------------------------
      if (! newcon) { ## 0 is FALSE and 1 is TRUE supposedly
         envjn$yk <- g - oldg 
         envjn$sk <- alpha*p 
         envjn$yksk <- as.numeric(crossprod(envjn$yk, envjn$sk)) 
         ## cat("yksk:",envjn$yksk,"\n")
         ireset <- ( (icycle == n-1) |  
                 (difnew < epsred*(fkeep-f)) ) 
         if (! ireset) { 
            envjn$yrsr <- 
                 as.numeric(crossprod(envjn$yr,envjn$sr)) 
            ireset <- (envjn$yrsr <= 0) 
         } 
         upd1 <- (envjn$yksk <= 0) 
      }
      ## cat("newcon, upd1:", newcon, upd1,"\n")
## ---------------------------------------------------------
##  compute the search direction
## ---------------------------------------------------------
      argvec <- c(accrcy, gnorm, xnorm) 
##   [p, gtp, ncg1, d] <- ...
##	   modlnp (d, x, g, maxit, upd1, ireset, bounds, 
##                   ipivot, argvec, sfun) 
      mres  <- modlnp (d, x, g, maxit, upd1, ireset,
                bounds, ipivot, argvec, sfun, ...) 
      ncg1 <- mres$ncg1
      gtp  <- mres$gtp
      d    <- mres$dnew
      p    <- mres$p

   ## cat("New p:")
   ## print(p)
   ## tmp<-readline("cont.")

      ncg <- ncg + ncg1 
## ---------------------------------------------------------
##  update LMQN preconditioner
## ---------------------------------------------------------
      if (! newcon) { 
         if (ireset) { 
            envjn$sr  <- envjn$sk 
            envjn$yr  <- envjn$yk 
            fkeep  <- f 
            icycle <- 1 
         } else {
            envjn$sr     <- envjn$sr + envjn$sk 
            envjn$yr     <- envjn$yr + envjn$yk 
            icycle <- icycle + 1 
         }
      }
   } # end while !conv
}
