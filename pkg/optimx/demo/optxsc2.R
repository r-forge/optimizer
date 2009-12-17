###################################################
## optxsc2 -- examples of use of optimx function
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

sc2.f <- function(x){
nvec <- 1:length(x)
sum(nvec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
nvec <- 1:length(x)
nvec * (exp(x) - 1) / 10
}
##-----------------------------------------------------##

#########################################################################################
p0 <- rnorm(20,sd=2)

ans.sc2df<-optimx(p0,fn=sc2.f,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem sc2 no gradients\n")
print(ans.sc2df)

ans.sc2g<-optimx(p0,fn=sc2.f,gr=sc2.g,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem sc2 with gradients\n")
print(ans.sc2g)

# This is to demonstrate the value of providing "exact" gradient information

##########
