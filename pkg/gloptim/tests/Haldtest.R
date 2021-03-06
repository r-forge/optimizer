##  Test function HALD, 5 parameters, non-smooth, one local minimum

# Hald <- list(
#     fn = function(x) {
#         stopifnot(is.numeric(x), length(x) == 5)
#         if (all(x == 1)) return(exp(1))
#         t <- -1 + (c(1:21) - 1)/10
#         v <- (x[1]+x[2]*t) / (1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
#         max(abs(v))
#     },
#     dim = 5, lb = rep(-1, 5), ub = rep(1, 5),
#     xmin = c(0.99987763, 0.25358844, -0.74660757, 0.24520150, -0.03749029), 
#     fmin = 0.00012237326, prec = 1.0e-10
# )
# 
# 
# mthd <- c(); fminimum <- c(); thetime = c()
# 
# for (m in c("deoptim", "cppdeoptim", "deoptimr",     # **DE**
#             "deopt", "simplede", "simpleea",         # **EA**
#             "gensa", "ga",                           # "gaopt" **GA**
#             "pso", "psopt", "hydropso",              # **PSO**
#             "ceimopt",                               # **CE**
#             "direct", "crs2lm", "isres",             # **NLoptr**
#             "cmaoptim", "cmaes", "purecmaes",        # **CMA-ES**
#             "malschains", "smco", "soma")) {
#     tm <- system.time(
#         sol <- gloptim(Hald$fn, rep(-1, 5), rep(1, 5), method=m)
#     )
#     mthd <- c(mthd, m); fminimum <- c(fminimum, sol$fmin)
#     thetime <- c(thetime, unname(tm["elapsed"]))
# }
# 
# G <- data.frame(method=mthd, fmin=fminimum, time=thetime)
# 
# G
# ##      method         fmin  time
# ##  1  deoptim 0.0001265921 1.034
# ##  2 deoptimr 0.0001223713 2.778
# ##  3 simplede 0.0001223713 2.909
# ##  4       ga 0.0117655364 5.275
# ##  5     smco 0.3015429636 0.053
# ##  6     soma 0.3772099183 0.116
# 
