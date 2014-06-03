## Run truncated-Newton method to optimize function represented 
## by grosesub.R

source("tn.R")
source("lmqn.R")
source("lin1.R")
source("modlnp.R")
source("gtims.R")
source("tn.R")
source("my.sfun.R")
source("initpc.R")
source("msolve.R")
source("ndia3.R")
source("step1.R")
source("ssbfgs.R")
source("grosesub.R")

## Initial guess

n <- 4 # normally 1000
x <- 5*c(pi/(1:n))


cat("Unconstrained example\n")
print(x)

ures  <- tn(x, grosefg)
print(ures)

