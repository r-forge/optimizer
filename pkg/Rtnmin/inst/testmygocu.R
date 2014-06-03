## Run truncated-Newton method to optimize function represented 
## by grosesub.R
require(optimx)
source("modz.R")
source("tn.R")
source("tnbc.R")
source("crash.R")
source("ztime.R")
source("stpmax.R")
source("cnvtst.R")
source("lmqn.R")
source("lmqnbc.R")
source("lin1.R")
source("modlnp.R")
source("gtims.R")
source("my.sfun.R")
source("initpc.R")
source("msolve.R")
source("ndia3.R")
source("step1.R")
source("ssbfgs.R")
source("grosesub.R")

## Initial guessrequire(optimx)
source("grosesub.R")

## Initial guess

n <- 4 # normally 1000
## x <- 5*rnorm(n)
x <- 5*c(pi/(1:n))
print(x)

cat("Unconstrained example \n")
cat("with tn\n")
ures <- tn(x, grosefg, trace=TRUE)
print(ures)
cat("with optimx\n")
ureso <- optimx(x, genrose.f, genrose.g, control=list(trace=2, all.methods=TRUE))
print(summary(ureso, order=value))

cat("Constrained example \n")
low <- -.02*rep(1,n) + .005*runif(n,1) ## lower bounds
up  <-  1.90*rep(1,n) + .005*runif(n,1) ## lower bounds

# to check -- run with very wide bounds

low<- -1e100*rep(1,n)
up<- -low

x<-0.5*(low+up)

print(low)
print(up)

cat("with tnbc\n")
cres1 <- tnbc(x, grosefg, low, up)
print(cres1)
cat("with optimx\n")
creso1 <- optimx(x, genrose.f, genrose.g, lower=low, upper=up, control=list(trace=2, all.methods=TRUE))
print(summary(creso1, order=value))


low<- 0:3
up <- 4:7

x<-0.5*(low+up)

print(low)
print(up)

cat("with tnbc\n")
cres2 <- tnbc(x, grosefg, low, up)
print(cres2)
cat("with optimx\n")
creso2 <- optimx(x, genrose.f, genrose.g, lower=low, upper=up, control=list(trace=2, all.methods=TRUE))
print(summary(creso1, order=value))
