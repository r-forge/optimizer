modlnp <- function(d, x, g, maxit, upd1, ireset, bounds, 
         ipivot, argvec, sfun, ...) {
##---------------------------------------------------------
## this routine performs a preconditioned conjugate-gradient
## iteration to solve the Newton equations for a search
## direction for a truncated-newton algorithm. 
## When the value of the quadratic model is sufficiently 
## reduced, the iteration is terminated.
##---------------------------------------------------------
## parameters
##
## p           - computed search direction
## g           - current gradient
## gv,gz1,v    - scratch vectors
## r           - residual
## d           - diagonal preconditoning matrix
## feval       - value of quadratic function
##------------------------------------------------------------
## initialization
##------------------------------------------------------------

   if (is.null(d)) stop("Null d")
   ## print(x)
   ## print(g)
   ## cat(bounds,"\n")

   accrcy <- argvec[1] 
   gnorm  <- argvec[2] 
   xnorm  <- argvec[3] 

   if (maxit == 0) {
      p    <- -g 
      gtp  <- as.numeric(crossprod(p, g))
      ncg1 <- 1
      dnew <- d
      if (sqrt(sum(p^2))==0) {
         ## cat("MODLNP 01: |p| = 0\n") 
         ## pause(1) 
      }
#%       cat("modlnp - dout01:")
#%       print(dnew)
      result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
      return(result) 
   }

   first <- 1
   tol   <- 1e-6   ######### was 1.d-12  #########
   qold  <- 0
   ainit  <- initpc (d, upd1, ireset)
   dnew<-ainit$td
   
#%    cat("after initpc dnew:")
#%    print(dnew)
   r     <- -g
   v     <- rep(0, length(r))
   p     <- v
   gtp   <- as.numeric(crossprod(p, g))

   rho  <- rep(0, maxit+1)
   beta <- rep(0, maxit)
   v.gv <- rep(0, maxit)

   rho[[1]] <- as.numeric(crossprod(r))

##------------------------------------------------------------
## main iteration (conjugate-gradient iterations for Ax = b)
##------------------------------------------------------------
   ind <- 0
   ncg1 <- 0
   for (k in 1:maxit) {
      ncg1 <- ncg1 + 1
      if (bounds) { r  <- ztime(r, ipivot) }
      amsolve <- msolve (r, upd1, ireset, first, d) 
      zk<-amsolve$y
      
      if (bounds) { zk <- ztime (zk, ipivot) }
#%       cat("r:")
#%       print(r)
#%       cat("zk:")
#%       print(zk)
      rz <- as.numeric(crossprod(r,zk))  

      if (rz/gnorm < tol) {
         ind <- 80 
         if (sqrt(sum(p^2))==0) {
            p <- -g
            gtp <- as.numeric(crossprod(p, g))
         }
#%          cat("modlnp - dout - ind80:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      } 
      if (k > 1) {
         beta[k] <- rz/rzold 
      } else {
         beta[k] <- 0 
      }
#%       cat("beta[",k,"] =",beta[k],"\n")
      v <- zk + beta[k]*v 
      if (bounds) { v  <- ztime( v, ipivot) }
      ## cat("about to call gtims\n")
      gv <- gtims(v, x, g, accrcy, xnorm, sfun, ...) 
#%     cat("After gtims: gv=")
#%     print(gv)
      ## cat(bounds,"\n")
      if (bounds) { gv <- ztime (gv, ipivot) }
#%       cat("gv, v:")
      ## cat("gv=", gv, "\n")
#%       print(gv)
#%       print(v)
      v.gv[[k]] <- as.numeric(crossprod(v, gv))  ## ?? v.gv has underscore!! ??
#%       cat("v.gv[[",k,"]]=",v.gv[[k]],"\n")
      if (v.gv[[k]]/gnorm < tol) { 
         ind <- 50  
         if (sqrt(sum(p^2))==0) {
            cat("MODLNP 03: |p| = 0 \n") 
            ## pause(1) 
         }
#%          cat("modlnp - dout03:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      }
      dnew <- ndia3(dnew, v, gv, r, v.gv[[k]]) # NB need to return something but it doesn't get used
#%       cat("after ndia3 below dout03, dnew:")
#%       print(dnew)
##------------------------------------------------------------
## compute current solution and related vectors
##------------------------------------------------------------
      alpha <- rz / v.gv[[k]]
      p <- p + alpha* v 
      r <- r - alpha*gv 

      rho[[k+1]] <- as.numeric(crossprod(r))

##------------------------------------------------------------
## test for convergence
##------------------------------------------------------------
      gtp <- as.numeric(crossprod(p,g))
      pr <-  as.numeric(crossprod(r,p))
      q <- (gtp + pr) / 2 
      qtest <- k * (1 - qold/q) 
      if (qtest <= 0.5)  {
         if (sqrt(sum(p^2))==0) {
            cat("MODLNP 04: |p| = 0\n") 
            ## pause(1) 
         }
#%          cat("modlnp - dout04:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      }
      qold <- q 
##------------------------------------------------------------
## perform cautionary test
##------------------------------------------------------------
      if (gtp > 0) {
         ind <- 40  
         if (sqrt(sum(p^2))==0) {
            cat("MODLNP 05: |p| = 0 \n") 
            ##   pause(1) 
         }
#%          cat("modlnp - dout05:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      } 
      rzold <- rz 
   } ## end loop over k
   k <- k-1 
##------------------------------------------------------------
## terminate algorithm
##------------------------------------------------------------
   if (ind == 40) {
      p   <- p - alpha*v 
   } 
##------------------------------------------------------------
   if (ind == 50 && k <= 1) {
      amsolve <- msolve (g, upd1, ireset, first, d) 
      p<-amsolve$y
      
      p <- -p 
      if (bounds) { p <- ztime (p, ipivot) }
   }
##------------------------------------------------------------
   if (ind == 80 && k <= 1) {
      p <- -g 
      if (bounds) { p <- ztime (p, ipivot) }
   }
##------------------------------------------------------------
## store new diagonal preconditioner
##------------------------------------------------------------
   gtp  <- as.numeric(crossprod(p, g)) 
   ncg1 <- k + 1 
   if (sqrt(sum(p^2))==0) { 
       cat("MODLNP 06: |p| = 0 \n")
       ##    pause(1) 
   }
   ## cat("modlnp - dout06:")
#%    print(dnew)
   result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
   return(result) 
}
