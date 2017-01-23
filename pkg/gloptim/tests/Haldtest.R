##  Test function HALD, 5 parameters, non-smooth, one local minimum
require(gloptim) # Is this necessary?

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


with(Hald, {
  stime = system.time(
    sol <- gloptim(fn = fn, lb = lb, ub = ub, method = "ga",
                   minimize = TRUE,
                   control = list(popsize = 200, itermax = 1000))
    )
  cat("xmin: ", sol$xmin, '\n')
  cat("fmin: ", sol$fmin, '\n')
  cat("xerr: ", sqrt(sum((sol$xmin-Hald$xmin)^2)), '\n')
  cat("ferr: ", abs( sol$fmin-Hald$fmin), '\n')
  cat("Elapsed time: ", stime["elapsed"], " [s].")
})
## Global solver/method: ga 
## xmin:  0.9747484 0.4539033 -0.5457178 0.02722956 0.04386974 
## fmin:  0.02726052 
## xerr:  0.3677573 
## ferr:  0.02713814 
## Elapsed time:  10.632  [s].
with(Hald, {
  stime = system.time(
    sol <- gloptim(fn = fn, lb = lb, ub = ub, method = "deoptim",
                   minimize = TRUE,
                   control = list(itermax = 1000, info = FALSE))
    )
  cat("xmin: ", sol$xmin, '\n')
  cat("fmin: ", sol$fmin, '\n')
  cat("xerr: ", sqrt(sum((sol$xmin-Hald$xmin)^2)), '\n')
  cat("ferr: ", abs( sol$fmin-Hald$fmin), '\n')
  cat("Elapsed time: ", stime["elapsed"], " [s].")
})
## Global solver/method: deoptim 
## xmin:  0.9998748 0.2534681 -0.746735 0.2453215 -0.03752767 
## fmin:  0.0001251754 
## xerr:  0.0002157009 
## ferr:  2.802174e-06 
## Elapsed time:  1.035  [s].

