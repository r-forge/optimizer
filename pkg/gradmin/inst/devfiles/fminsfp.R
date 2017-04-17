fminsfp<-function (f, b, st, fmin=NULL, maxiter = 1000, maximum = FALSE, 
    tol = 1e-07, rel.tol = tol, abs.tol = 1e-15, ...) { # Nash alg 17 in R
  control=list(trace=2)
  fun <- match.fun(f)
  if (maximum) 
      f <- function(x) -fun(x, ...)
  else f <- function(x) fun(x, ...)
#    phi <- 0.5 * (3 - sqrt(5))
  big <- .Machine$double.xmax
  offset <- 100.0
  ifn <- 0
  if (is.null(fmin)) { # fmin is S1
      ifn <- ifn + 1
      P <- f(b)
      if (control$trace > 0) cat("top: f(",b,")=",fmin,"\n")
  } else { P <- fmin }
  smult <- 1.5 # A1
  fmult <- -0.25 # A2?? put in control list
  while (abs(st) > abs.tol) { # main loop
    S1 <- P
    S0 <- -big
    x1 <- 0
    bmin <- b
    repeat { # loop until we have S0 >= S1 < P (a valley)
      x2 <- x1 + st # step 3
      b <- bmin + x2
      if ((offset + b) == (offset + bmin + x1)) {
        res <- bmin
        attr(res, "fmin") <- S1
        attr(res, "fcount") <- ifn
        cat("Finished with b = bmin + x1\n")
        return(res)
      }
      ifn <- ifn + 1
      P <- f(b)
      cat("s/f eval: f(",b,")=",P,"\n")
      if (P < S1) { # success
        cat("Success\n")
        x0 <- x1
        S0 <- S1
        x1 <- x2
        S1 <- P
        st <- smult * st # then go to step 3 and retry
      } else { # failure
        cat("failure\n")
        if (S0 >= S1) break # we have valley condition
        S0 <- P
        x0 <- x2
        st <- fmult * st
      }
      tmp <- readline("End s/f repeat")
    } # end repeat loop
    # Parabolic Inverse Interpolation
    cat("Paramin\n")
    x0 <- x0 - x1
    S0 <- (S0 - S1)*st
    P <- (P - S1)*x0
    if ((offset + P) == (offset + S0)) { # paramin failure case
       b <- bmin + x1
       P <- S1
    } else {
       st <- 0.5*(P*x0 - S0*st)/(P - S0)
       x2 <- x1 + st
       b <- bmin + x2
       if ((offset + b) != (offset + x1 + bmin)) {
          ifn <- ifn + 1
          P <- f(b)
          cat("Par point: f(",b,")=",P,"\n")
          if (P < S1) {
            x1 <- x2 # step to minimum position
          } else {
            b <- bmin + x1
            P <- S1
          }
       }
    }
    st <- fmult*st # to force new try
    tmp <- readline("End main loop")
  } # end main loop
  return(list(xmin = bmin, fmin = S1, niter = ifn, estim.prec = abs(st)))
}
