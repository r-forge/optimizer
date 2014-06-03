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

library(microbenchmark)

n <- 1000 # normally 1000
set.seed(123.456)
x <- abs(2+2*rnorm(n))

cat("Unconstrained example\n")
print(x)

tures<-microbenchmark(ures  <- tn(x, grosefg, trace=FALSE))
print(ures)
print(tures)
print(memory.profile())

library(Rcgmin)
tcg<-microbenchmark(cgres<-Rcgminu(x, genrose.f, genrose.g))
print(cgres)
print(tcg)


cat("mean time cgres =",mean(tcg$time)," to fn ",cgres$f,"\n")
cat("gnorm=",as.numeric(crossprod(genrose.g(cgres$par))),"\n")

library(nloptwrap)

nl.opts(list(maxeval=1000*n))
ttnnl<-microbenchmark(tnnlres<-tnewton(x, genrose.f, genrose.g))
print(tnnlres)
print(ttnnl)


cat("mean time tnnlres =",mean(ttnnl$time)," to fn ",tnnlres$value,"\n")
cat("gnorm=",as.numeric(crossprod(genrose.g(tnnlres$par))),"\n")

