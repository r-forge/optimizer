fuzz <- 1e-1 #1e-3 # 1e-10
options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("BB test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

p0 <- runif(500,min=-1, max=1)
system.time(ans.spg <- spg(par=p0, fn=sc2.f, control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
good <- -0.02404571778303681
print(z, digits=16)
if(any(abs(good - z) > fuzz)) stop("BB test sc2 a FAILED")

system.time(ans.spg <- spg(par=p0, fn=sc2.f, grad=sc2.g,
   control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
good <- 0.000849270257616394
print(z, digits=16)
if(any(abs(good - z) > fuzz)) stop("BB test sc2 b FAILED")

system.time(ans.opt <- optim(par=p0, fn=sc2.f, method="L-BFGS-B"))[1]

z <- sum(ans.opt$par)
good <- 0.02209066162550582
print(z, digits=16)
#if(any(abs(good - z) > fuzz)) stop("BB test sc2 c FAILED")

system.time(ans.opt <- optim(par=p0, fn=sc2.f, gr=sc2.g, method="L-BFGS-B"))[1]

z <- sum(ans.opt$par)
good <- 0.02200130759852783
print(z, digits=16)
if(any(abs(good - z) > fuzz)) stop("BB test sc2 d FAILED")

