cgsol24 <- function(AA, c) {
  # Algorithm 24 conjugate gradient solution of equations for posdef symm AA
  # Need a tol scaled to problem
  n <- length(c)
  b <- rep(0, n)
  tol <- .Machine$double.eps * (max(AA)+max(b))/2
  offset <- 1e+6
  cat("Tolerance to use in cgsol24=", tol, "\n")
  # Step 1
#  g <- ax(b, AA) - c
  g <- c
  G2 <- as.numeric(crossprod(g))
  t <- -g
  for (itn in 1:n){
     cat("itn ",itn,"  G2 =",G2,"\n")     
     if (G2 <= tol) break
     v <- ax(t,AA)
     T2 <- as.numeric(crossprod(t, v))
     k <- G2/T2
     G2L <- G2
     g <- g + k*v
     G2 <- as.numeric(crossprod(g))
     bt <- b + k*t
     if (sum(((bt+offset)==(b+offset))) == n) break
     b <- bt
     T2 <- G2/G2L
     t <- T2*t - g
  }
  # Do we want to restart?
  list(x=b, ressq=G2)
}