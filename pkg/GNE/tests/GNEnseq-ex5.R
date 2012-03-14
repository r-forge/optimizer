
require(GNE)	



#-------------------------------------------------------------------------------
# (4) Example of GNE with 1 solution
#-------------------------------------------------------------------------------



F <- function(z, phi=phiMin, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	# cat("x", x, "\n")
	c(	-2*(x[1] - 2)/x[2] + lambda[1],
		-2*(x[2] - 3)/x[1] + lambda[2],
		phi(lambda[1], x[1], ...),
		phi(lambda[2], x[2], ...)
	)
}


JacF <- function(z, gphia, gphib, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	idga1 <- gphia(lambda[1], x[1], ...)
	idgb1 <- gphib(lambda[1], x[1], ...)
	idga2 <- gphia(lambda[2], x[2], ...)
	idgb2 <- gphib(lambda[2], x[2], ...)

rbind(
c(-2, 2*(x[1] - 2)/x[2]^2, 1, 0),
c(2*(x[2] - 3)/x[1]^2, -2, 0, 1),
c(idgb1, idgb1, idga1, 0),
c(idgb2, idgb2, 0, idga2)
)	
}



#a simple test
z0 <- c(3, 4, 1, 1)

F(z0)
JacF(z0, gphia= GrAphiFB, gphib= GrBphiFB)

GNE.nseq(z0, F, JacF, list(phi=phiMin), list(gphia= GrAphiMin, gphib= GrBphiMin), method="Newton")

GNE.nseq(z0, F, JacF, list(phi= phiFB), list(gphia= GrAphiFB, gphib= GrBphiFB), method="Newton")

#-------------------------------------------------------------------------------
#random initial points
n <- 20
set.seed(1234)
initpt <- cbind(runif(n, -10, 10), runif(n, -10, 10), 1, 1)

NewLnsrch <- function(i, echo=TRUE, phi= phiFB, gphia= GrAphiFB, gphib= GrBphiFB)
{ 
	if(echo)
		cat("______", initpt[i, ], "\n")
	res <- GNE.nseq(initpt[i, ], F, JacF, list(phi=phi), list(gphia=gphia, gphib=gphib), method="Newton") 
	if(echo)
		print(res$par)
	c(res$par, res$value )
}

respt <- t( sapply(1:NROW(initpt), NewLnsrch ) )

totalres <- cbind(1:NROW(initpt), initpt[, 1:2], NA, round(respt, 3))

#remove non convergent optimization sequences
finalres <- totalres[totalres[,9] == 0, -9]
colnames(finalres) <- c("num", "x1 init", "x2 init", "", "x1 final", "x2 final", "lam 1", "lam 2")

print(finalres)




	