###################################################
## optxrosbkext -- examples of use of optimx function
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

rosbkext.f <- function(x){
p <- x
n <- length(p)
sum (100*(p[1:(n-1)]^2 - p[2:n])^2 + (p[1:(n-1)] - 1)^2)
}
##-----------------------------------------------------##


##########
p0 <- rnorm(20,sd=2)
ans.rosbext<-optimx(p0,fn=rosbkext.f,control=list(all.methods=TRUE, save.failures=TRUE, maxit=500))
cat("\n")
cat("Problem rosbext\n")
print(ans.rosbext)
 
