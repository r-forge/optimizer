## Run truncated-Newton method to optimize function represented 
## by grosesub.R

source("tnbc.R")
source("lmqnbc.R")
source("crash.R")
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
source("ztime.R")
source("stpmax.R")
source("cnvtst.R")
source("modz.R")
source("norm2.R")

## Initial guess

n <- 4 # normally 1000
## x <- 5*rnorm(n)
x <- 5*c(pi/(1:n))
print(x)

cat("Constrained example \n")
low <- -.02*rep(1,n) + .005*runif(n,1) ## lower bounds
up  <-  1.90*rep(1,n) + .005*runif(n,1) ## lower bounds

# to check -- run with very wide bounds

low<- -1e100*rep(1,n)
up<- -low

low<- 0:3
up <- 4:7

x<-0.5*(low+up)


print(low)
print(up)
print(x)

cres <- tnbc(x, grosefg, low, up)
print(cres)
