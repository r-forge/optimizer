



partsvd <- function(A, ns=25, cyclelimit=6) {
# Partial svd by the one-sided Jacobi method of Nash & Shlien
#  Computer Journal 1987 ???
  eps <- .Machine$double.eps
  if (cyclelimit < 6) cyclelimit <- 6 # safety in case user tries smaller
  m <- dim(A)[1]
  n <- dim(A)[2]
  e2 <- 10*m*eps*eps
  tol <- 0.1 * eps
  EstColRank <- n # estimated column rank 
  # Note that we may simply run algorithm to completion, or fix the
  # number of columns by ns. Need ?? to fix ns=0 case.??
  V <- diag(nrow=n) # identity matrix in V
  if (is.null(ns)) {ns <- n } # Safety check on number of svs
  z <- rep(NA, n) # column norm squares -- safety setting
  keepgoing <- TRUE
  SweepCount <- 0
  while (keepgoing) { # main loop of repeating cycles of Jacobi
    RotCount <- EstColRank*(EstColRank - 1)/2
    SweepCount <- SweepCount + 1 
    ns <- EstColRank
    for (jj in 1:(ns - 1)) { # left column indicator
       for (kk in (jj+1): ns) { # right hand column
         p <- q <- r <- 0.0 # 
         p <- as.numeric(crossprod(A[,jj], A[,kk]))
         q <- as.numeric(crossprod(A[,jj], A[,jj]))
         r <- as.numeric(crossprod(A[,kk], A[,kk]))
         z[jj]<-q
         z[kk]<-r
         if (q >= r) { # in order, so can do test of "convergence"
            if ( (q <= e2*z[1]) || (abs(p) <= tol*q) ) {
                RotCount <- RotCount - 1 # ignore rotation
                break # hopefully goes to next kk
            }
            p <- p/q
            r <- 1 - (r/q)
            vt <- sqrt(4*p*p +r*r)
            c0 <- sqrt(0.5*(1+r/vt))
            s0 <- p/(vt*c0)
            # rotate
            cj <- A[,jj]
            ck <- A[,kk]
            A[,jj] <- c0*cj + s0*ck
            A[,kk] <- -s0*cj + c0*ck
            cj <- V[,jj]
            ck <- V[,kk]
            V[,jj] <- c0*cj + s0*ck
            V[,kk] <- -s0*cj + c0*ck
         } else { # out of order, must rotate
            p <- p/r
            q <- (q/r) - 1.0
            vt <- sqrt(4*p*p +q*q)
            s0 <- sqrt(0.5*(1-q/vt))
            if (p < 0) { s0 <- -s0 }
            c0 <- p/(vt*s0)
            # rotate
            cj <- A[,jj]
            ck <- A[,kk]
            A[,jj] <- c0*cj + s0*ck
            A[,kk] <- -s0*cj + c0*ck
            cj <- V[,jj]
            ck <- V[,kk]
            V[,jj] <- c0*cj + s0*ck
            V[,kk] <- -s0*cj + c0*ck
         } # end q >= r test
       } # end kk
    } # end jj
    cat("End sweep ", SweepCount,"  No. rotations =",RotCount,"\n")
    while( (EstColRank >= 3) && (z[EstColRank] <= (z[1]*tol+tol*tol)) ) {
    # ?? Why can we not use 2?
        EstColRank <- EstColRank - 1
        cat("Reducing rank to ", EstColRank,"\n") # ?? can do this more cleanly
    } # end while for rank estimation
    if ( SweepCount >= cyclelimit) { 
         cat("Cycle limit reached\n")
         keepgoing <- FALSE
    } 
    if (RotCount == 0) {
        cat("Zero rotations in cycle\n")
        keepgoing <- FALSE
    }
  } # End main cycle loop
  ans <- list( sqsv = z, U = A, V=V)
  ans
} # end partsvd()
 
