##  Test function HALD, 5 parameters, non-smooth, one local minimum
## Haldtestj.R -- for JN functions

# Hald <- list(
#   fn = function(x) {
#     stopifnot(is.numeric(x), length(x) == 5)
#     if (all(x == 1)) return(exp(1))
#     t <- -1 + (c(1:21) - 1)/10
#     v <- (x[1]+x[2]*t) / (1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
#     max(abs(v))
#   },
#   dim = 5, lb = rep(-1, 5), ub = rep(1, 5),
#   xmin = c(0.99987763, 0.25358844, -0.74660757, 0.24520150, -0.03749029), 
#   fmin = 0.00012237326, prec = 1.0e-10
# )
# 
# 
# with(Hald, {
#   stime = system.time(
#     sol <- gloptimj(fn = fn, lb = lb, ub = ub, method = "smco",
#                    minimize = TRUE,
#                    control = list(itermax = 1000, info = FALSE))
#   )
#   cat("xmin: ", sol$xmin, '\n')
#   cat("fmin: ", sol$fmin, '\n')
#   cat("xerr: ", sqrt(sum((sol$xmin-Hald$xmin)^2)), '\n')
#   cat("ferr: ", abs( sol$fmin-Hald$fmin), '\n')
#   cat("Elapsed time: ", stime["elapsed"], " [s].")
# })
# ## Global solver/method: smco 
# ## xmin:  0.9998748 0.2536428 -0.7465512 0.2451177 -0.0374414 
# ## fmin:  0.000126864 
# ## xerr:  0.0001246979 
# ## ferr:  4.490704e-06 
# ## Elapsed time:  3.361  [s].  
# 
# with(Hald, {
#   stime = system.time(
#     sol <- gloptimj(fn = fn, lb = lb, ub = ub, method = "soma",
#                    minimize = TRUE,
#                    control = list(itermax = 1000, info = FALSE))
#   )
#   cat("xmin: ", sol$xmin, '\n')
#   cat("fmin: ", sol$fmin, '\n')
#   cat("xerr: ", sqrt(sum((sol$xmin-Hald$xmin)^2)), '\n')
#   cat("ferr: ", abs( sol$fmin-Hald$fmin), '\n')
#   cat("Elapsed time: ", stime["elapsed"], " [s].")
# })
# ## Global solver/method: soma 

