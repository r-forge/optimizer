fuzz <- 1e-3 # 1e-10
options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test chen.f ...\n")

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

p0 <- rexp(500)
system.time(ans.spg <- spg(par=p0, fn=chen.f, lower=0))[1]
system.time(ans.opt <- optim(par=p0, fn=chen.f, lower=0, method="L-BFGS-B"))[1]
 
z <- sum(ans.spg$par)
good <- 1801.831961319466
print(z, digits=16)
if(any(abs(good - z)/1e4 > fuzz)) stop("BB test chen.f a FAILED")
 
z <- sum(ans.opt$par)
good <- 2243.132018091285
print(z, digits=16)
if(any(abs(good - z)/1e4 > fuzz)) stop("BB test chen.f b FAILED")

