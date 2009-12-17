###################################################
## optxtrig -- examples of use of optimx function
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


trig.f <- function(x){
n <- length(x)
i <- 1:n
f <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
sum(f*f)
}
##-----------------------------------------------------##

 
##########
p0 <- rnorm(20,sd=5)
ans.trig<-optimx(p0,fn=trig.f,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem trig\n")
print(ans.trig)

