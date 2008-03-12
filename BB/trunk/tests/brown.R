options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test brown.f ...\n")

brown.f <- function(x) {
p <- x
n <- length(p)
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}

p0 <- rnorm(500,sd=2) # this set fails in optim, so
p0 <- rnorm(500,sd=2)
system.time(ans.spg <- spg(par=p0, fn=brown.f, maxit=2500))[1]

z <- sum(ans.spg $par)
good <- -0.02486011469681782
print(z, digits=16)
if(any(abs(good - z) > 1e-10)) stop("BB test brown.f a FAILED")

system.time(ans.opt <- optim(par=p0, fn=brown.f, method="L-BFGS-B"))[1]

z <- sum(ans.opt $par)
good <- 0.00898911426823444
print(z, digits=16)
if(any(abs(good - z) > 1e-10)) stop("BB test brown.f b FAILED")

