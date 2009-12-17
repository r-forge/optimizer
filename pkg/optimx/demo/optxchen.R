###################################################
## optxchen -- examples of use of optimx function
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

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}
##-----------------------------------------------------##


##########
p0 <- rexp(50)
ans.chen<-optimx(p0,fn=chen.f,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem chen\n")
print(ans.chen)

