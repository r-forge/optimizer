###################################################
## optxmin -- examples of use of optimx function
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


froth <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
sum (f * f)
}
##-----------------------------------------------------##


#########################################
p0 <- rpois(2,10)
ans.froth<-optimx(p0,fn=froth,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem froth\n")
print(ans.froth)


