if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1234)
old.seed <- setRNG(test.rng)


# Nonlinear resistive circuit (Yamamura, JCAM 2003)
nrc <- function(x) {
  	n <- length(x)
  	f <- rep(NA, n)
	f <- 2.5*x^3 - 10.5*x^2 + 11.8*x + sum(x) - (1:n)
	f
}

p0 <- matrix(rnorm(500), 50, 10)  # 100 starting values, each of length 5
#ans1 <- dfsane(par=p0[1,], fn=nrc)
#ans <- BBsolve(par=p0, fn=nrc)
ans <- multiStart(par=p0, fn=nrc)

#pc <- princomp(ans$par) 
#plot(pc$scores[,1])  # plot shows distinct solutions


# ans$convergence
 
z <- sum(ans$par)
good   <-    265.9957598676737  # early defaults: 267.0
#on Windows 
#on Linux64  early defaults: 266.8965710381518
#on Linux32  early defaults: 267.2236080079449
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test BBsolve NRC FAILED")
