## Run truncated-Newton method to optimize function represented 
## by my_sfun.R

source("tn.R")
source("tnbc.R")
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
source("lmqnbc.R")
source("crash.R")

## Initial guess

n <- 3 # normally 1000
x <- 1:n # was rnorm(n)

cat("Unconstrained example\n")
print(x)

ures  <- tn(x, my.sfun)
print(ures)

 cat("Constrained example \n")

## low <- -.02*rep(1,n) + .005*randn(n,1) ## lower bounds
## up  <-  .90*rep(1,n) + .005*randn(n,1) ## lower bounds
 low <- -.02*rep(1,n)  ## lower bounds
 up  <-  .90*rep(1,n)  ## lower bounds

 cres <- tnbc(x, my.sfun, low, up)
 print(cres)
