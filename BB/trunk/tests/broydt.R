options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test broydt.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

p0 <- rnorm(100, sd=1)
system.time(ans.spg <- spg(p=p0, func=broydt.f))[1]
system.time(ans.opt <- optim(par=p0, fn=broydt.f, method="L-BFGS-B"))[1]
 
z <- sum(ans.spg$par)
good <-  137.6240252276294
print(z, digits=16)
if(any(abs(good - z) > 1e-10)) stop("BB test broydt.f a FAILED")
 
z <- sum(ans.opt$par)
good <- 111.5078705487698
print(z, digits=16)
if(any(abs(good - z) > 1e-10)) stop("BB test broydt.f b FAILED")
