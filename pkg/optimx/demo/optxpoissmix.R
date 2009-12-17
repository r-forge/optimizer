###################################################
## optxpoissmix -- examples of use of optimx function
library(optimx)

#NOTE: bobyqa is failing!!

if(!require("BB"))    stop("this requires package BB.")
if(!require(numDeriv))stop("this requires package numDeriv.")
if(!require("setRNG"))stop("this requires setRNG.")
if(!require("Rcgmin"))stop("this requires Rcgmin.")
if(!require("minqa"))stop("this requires minqa.")
if(!require("ucminf"))stop("this requires ucminf.")

# Replaced April 7, 2008, with setRNG to ensure rng and normal generators are set too.
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=1236)
setRNG(test.rng)

###############################################################
poissmix.loglik <- function(p,y) {
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return (sum(loglik) ) # temp fix
}

###############################################################
# Real data from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

lo <- c(0.001,0,0)
hi <- c(0.999, Inf, Inf)

y <- poissmix.dat$freq

p0 <- runif(3,c(0.2,1,1),c(0.8,5,8))  # randomly generated starting values

cat("p0\n")
print(p0)

cat("check initial fn value: ")
fninit<-poissmix.loglik(p0,y=y)
cat(fninit,"\n")
#anspoiss.nm<-optimx(p0, fn=poissmix.loglik, y=y, method='Nelder-Mead', control=list(trace=2, fnscale=-1))
#tmp<-readline("OK")

ans.poiss1<-optimx(p0, fn=poissmix.loglik, y=y, lower=lo, upper=hi, control=list(maximize=TRUE, all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem poissmix\n")
print(ans.poiss1)

