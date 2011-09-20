library(GNE)

#constants
d <- 20
lambda <- 4
rho <- 1

#Phi(z) function  
phiex2 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]	
	c(	-(d - lambda - rho*(x+y) - rho*x) - s1,
		-(d - lambda - rho*(x+y) - rho*y) - s2,
		min(x, s1),
		min(y, s2)	)
}

#Jac Phi(z) function
jacphiex2 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]

	res <- matrix(0, 4, 4)
	res[1, ] <- c(2*rho, rho, -1, 0)
	res[2, ] <- c(rho, 2*rho, 0, -1)
	res[3, ] <- c( 1*(x <= s1), 0, 1*(s1 <= x), 0)
	res[4, ] <- c( 0, 1*(y <= s2), 0, 1*(s2 <= y))
	res
	
}

#call, true value is (16/3, 16/3) 

NewtonKKT(rep(0, 4), "Levenberg-Marquardt", phiex2, jacphiex2, control=list(maxit=200, echo=TRUE))

NewtonKKT(rep(0, 4), "Levenberg-Marquardt", phiex2, jacphiex2, control=list(maxit=200, echo=2))

NewtonKKT(rep(0, 4), "Levenberg-Marquardt", phiex2, jacphiex2, control=list(maxit=200, echo=3))


NewtonKKT(rep(0, 4), "New", phiex2, jacphiex2, control=list(maxit=200, echo=TRUE))


NewtonKKT(rep(0, 4), "Leven", phiex2, jacphiex2, control=list(maxit=200, echo=TRUE), extrapol="MPE")


NewtonKKT(rep(0, 4), "Leven", phiex2, jacphiex2, control=list(maxit=200, echo=TRUE), extrapol="SqMPE")



