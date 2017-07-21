# genrosa function code -- attempts to match the rosenbrock at gs=100 and x=c(-1.2,1)
genrosa.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
        # Note do not at 1.0 so min at 0
	fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

genrosa.g <- function(x, gs=NULL){
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

genrosa.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
#		z2<-1.0 - x[i-1]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

require(snewton)
cat("Generalized Rosenbrock tests\n")

cat("original function")

x0 <- c(-1.2, 1)
solorig <- snewton(x0, genrosa.f, genrosa.g, genrosa.h)
print(solorig)
print(eigen(solorig$Hess)$values)
solorigm <- snewtonm(x0, genrosa.f, genrosa.g, genrosa.h)
print(solorigm)
print(eigen(solorigm$Hess)$values)

cat("Start with 50 values of pi and scale factor 10\n")
x0 <- rep(pi, 50)
sol50pi <- snewton(x0, genrosa.f, genrosa.g, genrosa.h, gs=10)
print(sol50pi)
print(eigen(sol50pi$Hess)$values)
sol50pim <- snewtonm(x0, genrosa.f, genrosa.g, genrosa.h, gs=10)
print(sol50pim)
print(eigen(sol50pim$Hess)$values)
