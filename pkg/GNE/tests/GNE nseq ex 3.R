
library(GNE)	



#-------------------------------------------------------------------------------
# (3) Example 5.3 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#constants
cstC <- cbind(c(.1, .12, .15), c(.01, .05, .01))
cstU <- cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75))
cstK <- c(100, 100)
cstE <- c(.5, .25, .75)
cstD <- c(3, .01)

constr <- function(x, index)
	sum(cstU[, index] * cstE * x) - cstK[index]

gradL <- function(z, index)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	
	p1 <- -( cstD[1] - cstD[2]*sum(x) - cstC[index, 1] - 2*cstC[index, 2]*x[index] - cstD[2]*x[index])
	p2 <-  mu1*cstU[index,1]*cstE[index] + mu2*cstU[index,2]*cstE[index] 
	p1 + p2
}	

#phi function
phiex3 <- function(z)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	
c(	gradL(z, 1),
	gradL(z, 2),
	gradL(z, 3),
	min(-constr(x, 1), mu1),
	min(-constr(x, 2), mu2) )
}	

#Jacobian of the phi function
jacphiex3 <- function(z)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	A <- -constr(x, 1) <= mu1
	B <- -constr(x, 2) <= mu2
	
	res <- matrix(0, 5, 5)
	res[1, ] <- c(2*cstD[2]+2*cstC[1,2], cstD[2], cstD[2], cstU[1,1]*cstE[1], cstU[1,2]*cstE[1])
	res[2, ] <- c(cstD[2], 2*cstD[2]+2*cstC[2,2], cstD[2], cstU[2,1]*cstE[2], cstU[2,2]*cstE[2])
	res[3, ] <- c(cstD[2], cstD[2], 2*cstD[2]+2*cstC[3,2], cstU[3,1]*cstE[3], cstU[3,2]*cstE[3])
	res[4, ] <- c(-cstU[1,1]*cstE[1]*1*A, -cstU[2,1]*cstE[2]*1*A, -cstU[3,1]*cstE[3]*1*A, 1*(!A), 0)
	res[5, ] <- c(-cstU[1,2]*cstE[1]*1*B, -cstU[2,2]*cstE[2]*1*B, -cstU[3,2]*cstE[3]*1*B, 0, 1*(!B))

	res
}

#call, true value around (21.146, 16.027, 2.724, 0.574, 0.000)

z0 <- rep(0, 5)
GNE.nseq(z0, phiex3, jacphiex3, method="Newton")

GNE.nseq(z0, phiex3, jacphiex3, method="Broyden")



	