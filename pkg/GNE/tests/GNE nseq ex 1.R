
library(GNE)	


#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#Phi(z) function
phiex1 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]
	mu <- z[5]

	c(	2*(x-1) - s1 + mu,
		2*(y-1/2) - s2 + mu,
		min(x, s1),
		min(y, s2),
		min(1 - x - y, mu) )
}

#Jacobian of Phi(z)
jacphiex1 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]
	mu <- z[5]
	
	res <- matrix(0, 5, 5)
	res[1, ] <- c(2, 0, -1, 0, 1)
	res[2, ] <- c(0, 2, 0, -1, 1)
	res[3, ] <- c( 1*(x <= s1), 0, 1*(s1 <= x), 0, 0)
	res[4, ] <- c( 0, 1*(y <= s2), 0, 1*(s2 <= y), 0)
	res[5, ] <- c( -1*(1 - x - y <= mu), -1*(1 - x - y <= mu), 0, 0, -1*(mu <= 1 - x - y))
	res
}

#true value is (3/4, 1/4, 0, 0, 1/2)

z0 <- rep(0, 5)
GNE.nseq(z0, phiex1, jacphiex1, method="Newton")

GNE.nseq(z0, phiex1, jacphiex1, method="Broyden")



	