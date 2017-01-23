##  Test function HALD, 5 parameters, non-smooth, one local minimum
Hald <- list(
  fn = function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    if (all(x == 1)) return(exp(1))
    t <- -1 + (c(1:21) - 1)/10
    v <- (x[1]+x[2]*t) / (1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
    max(abs(v))
  },
  dim = 5, lb = rep(-1, 5), ub = rep(1, 5),
  xmin = c(0.99987763, 0.25358844, -0.74660757, 0.24520150, -0.03749029), 
  fmin = 0.00012237326, prec = 1.0e-10
)
