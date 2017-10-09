## @knitr BrownAlmostLinear.prb
# This is file BrownAlmostLinear.prb
probname <- "BrownAlmostLinear"
probdesc <- "This problem has variable number of parameters.
It is a nonlinear least squares problem with m = n, and
seems to be suitable for nonlinear equations, as there is
a solution of 0 at (alpha, alpha, ...., alpha, alpha^(1-n))
where alpha satisfies
n alpha^n - (n+1) alpha^(n-1) + 1 = 0 and in particular 
	alpha = 1.
f = 1 at (0,...,0, (n+1))
x0 = rep(0.5, n)
"
# Reference: BROWN, K.M. A quadratically convergent Newton-like 
# method based upon Gaussian elimination. 
# S I A M J. Numer. Anal. 6 (1969), 560-569.
# Also More et al. problem #27

balf.f <- function(x) {   ## Brown Almost Linear Function
	n <- length(x)
	f<-as.vector(matrix(0,n,1))
	sumx<-sum(x)
	for (i in 1:(n-1)) {
		f[i] <- x[i] + sumx - (n+1);
	}
	f[n] <- prod(x) -1;
  sum(f^2);
}

balf.g <- function(x) {
  n <- length(x)
  res <- balf.res(x)
  JJ <- balf.jac(x)
  gg <- 2*as.vector(t(JJ) %*% res)
  gg
}

balf.res <- function(x) {   ## Brown Almost Linear Function
	n <- length(x)
	f <- as.vector(matrix(0,n,1))
	sumx <- sum(x)
	for (i in 1:(n-1)) {
		f[i] <- x[i] + sumx - (n+1)
	}
	f[n] <- prod(x) - 1 
	as.vector(f)
}

balf.jac <- function(x){    ## Brown Almost Linear Function
  n <- length(x)
  J <- matrix(0, n, n)
  for (i in 1:(n-1)) {
     for (j in 1:n) {
        J[i,j] <- 1 # from sum(x)
     }
     J[i,i] <- 2 # 1 from x[i], 1 from sum(x)    
  }
  for (j in 1:n) {
     J[n,j] <- prod(x[-j])
  }
  attr(J,"gradient")<-J
  J
}

cat("Brown Almost Linear Problem\n")

n <- 10
start0 <- rep(0.5, n)
nlist <- rep("b",n)
for (i in 1:n) {
  nlist[i] <- paste(nlist[i],i,sep='')
}
names(start0) <- nlist
cat("start0:")
print(start0)

library(numDeriv)
x<-start0
cat("f:", balf.f(x),"\n")
cat("g:")
print(balf.g(x))
gn <- grad(balf.f, x)
cat("numgrad=")
print(gn)
cat("maxdiff =", max(abs(gn-balf.g(x))),"\n")

cat("res:")
rr <- balf.res(x)
print(rr)
cat("Jac:")
JA <- (balf.jac(x))
print(JA)
JN <- jacobian(balf.res, start0)
print (JN)
cat("maxdiff =", max(abs(JN-JA)),"\n")

fsol <- function(alpha, n) { 
  val <- n * alpha^n - (n+1) * alpha^(n-1) + 1 
  val
}

alstar <- uniroot(fsol, c(.8,.9999), n=n)$root
cat("Solution via uniroot for n=",n," has p1-(n-1) =", alstar," pn =", alstar^(1-n), "\n")

library(optimr)
BALP10opm <- opm(start0, balf.f, balf.g, method="ALL", control=list(kkt=FALSE))
summary(BALP10opm, order=value)

library(nlsr)

BALP10nlfb <- nlfb(start0, balf.res, jacfn=balf.jac, trace=TRUE)
print(BALP10nlfb)


