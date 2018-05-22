rm(list=ls())
require(optimx)
# genrose function code 
# A different generalization of the Rosenbrock banana valley function 
# -- attempts to match the rosenbrock at gs=100 and x=c(-1.2,1)
genrose.f<- function(x, gs=NULL){ # objective function
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
        # Note do not at 1.0 so min at 0
	fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn1]
        # f = gs*z1*z1 + z2*z2
	gg[tn] <- 2 * (gs * z1)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 - 2 *z2 
	return(gg)
}
xx <- rep(pi, 4)
lo <- rep(2, 4)
up <- rep(3, 4)

# !! NOTE: Parameter 4 bigger than upper bound
test1 <- optimr(xx, genrose.f, genrose.g, lower=lo, upper=up, method="Rvmmin", control=list(trace=3))
print(test1)
pstar <- test1$par
print(pstar)
pbetter<-pstar
pbetter[4]<-3
print(pbetter)
print(genrose.f(pstar))
print(genrose.f(pbetter))
test2 <- optimr(xx, genrose.f, genrose.g, lower=lo, upper=up, method="nlminb", control=list(trace=3))
print(test2)
test3 <- optimr(xx, genrose.f, genrose.g, lower=lo, upper=up, method="Rcgmin", control=list(trace=3))
print(test3)
test4 <- optimr(xx, genrose.f, genrose.g, lower=lo, upper=up, method="Rtnmin", control=list(trace=3))
print(test4)

require(Rvmmin)
test1a<-Rvmmin(xx, genrose.f, genrose.g, lower=lo, upper=up, control=list(trace=3))
print(test1a)
bc<-bmchk(xx, lower=lo, upper=up, shift2bound=TRUE)
test1b<-Rvmminb(bc$bvec, genrose.f, genrose.g, lower=lo, upper=up, bdmsk=bc$bdmsk, control=list(trace=3))
print(test1b)
# Note how you MUST start with parameters that are feasible!
