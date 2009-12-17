###################################################
## optxvalley -- examples of use of optimx function
library(optimx)

if(!require("BB"))    stop("this requires package BB.")
if(!require(numDeriv))stop("this requires package numDeriv.")
if(!require("setRNG"))stop("this requires setRNG.")
if(!require("Rcgmin"))stop("this requires Rcgmin.")
if(!require("minqa"))stop("this requires minqa.")
if(!require("ucminf"))stop("this requires ucminf.")

# Replaced April 7, 2008, with setRNG to ensure rng and normal generators are set too.
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=1236)
setRNG(test.rng)

## Functions for the tests **


valley.f <- function(x) {
c1 <- 1.003344481605351
c2 <- -3.344481605351171e-03
n <- length(x)
f <- rep(NA, n)
j <- 3 * (1:(n/3))
jm2 <- j - 2
jm1 <- j - 1
f[jm2] <- (c2 * x[jm2]^3 + c1 * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
f[j] <- 10 * (cos(x[jm2]) - x[j])
sum(f*f)
}
##-----------------------------------------------------##
 
##########
p0 <- rnorm(99, sd=2)
ans.valley<-optimx(p0,fn=valley.f,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem valley\n")
print(ans.valley)

