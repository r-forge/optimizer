
library(GNE)


#---------------------------------------------------------------
# example with 4 GNEs
#---------------------------------------------------------------

F <- function(z)
{
	x <- z[1:2]
	lambda <- z[3:4]
	# cat("x", x, "\n")
	c(	2*(x[1] - 2)*(x[2]-4)^4 + lambda[1],
		2*(x[2] - 3)*x[1]^4 + lambda[2],
		min(lambda[1], 1-sum(x)),
		min(lambda[2], 2-2*x[1]-x[2])
	)
}


JacF <- function(z)
{
	x <- z[1:2]
	lambda <- z[3:4]
	idg1 <- 1-sum(x) <= lambda[1]
	idg2 <- 2-2*x[1]-x[2] <= lambda[2]
rbind(
c(2*(x[2]-4)^4, 8*(x[1] - 2)*(x[2]-4)^3, 1, 0),
c(8*(x[2] - 3)*x[1]^3, 2*x[1]^4, 0, 1),
c(-1*idg1, -1*idg1, 1*!idg1, 0),
c(-2*idg2, -1*idg2, 0, 1*!idg2)
)	
}


args(GNE)

z0 <- c(10, Inf, 1, 1)



GNE("non smooth", "default", z0, F, JacF)

z0 <- c(10, 10, 1, 1)

GNE("non smooth", "default", z0, F, JacF)

GNE("non smooth", "default", z0, F, JacF, global="qline")

GNE("non smooth", "Broyden", z0, F, JacF, global="qline")


GNE.nseq(z0, F, JacF)
GNE.nseq(z0, F, JacF, method="Newton", global="dbldog")


?GNE
?GNE.nseq

