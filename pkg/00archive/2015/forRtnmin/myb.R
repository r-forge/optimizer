## Run truncated-Newton method to optimize function represented 
## by grosesub.R

source("crash.R")
source("lmqnbc.R")
source("lin1.R")
source("modlnp.R")
source("gtims.R")
source("tnbc.R")
source("my.sfun.R")
source("initpc.R")
source("msolve.R")
source("ndia3.R")
source("step1.R")
source("stpmax.R")
source("ssbfgs.R")
source("grosesub.R")

## Initial guess

n <- 4 # normally 1000
x <- 5*rnorm(n)

cat("Constrained example \n")

low <- -.02*rep(1, n) + .005*runif(n) ## lower bounds
up  <-  .90*rep(1, n) + .005*runif(n) ## upper bounds

cres <- tnbc(x, my.sfun, low, up)
print(cres)
