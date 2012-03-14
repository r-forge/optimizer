
require(GNE)	



#-------------------------------------------------------------------------------
# (4) Example of GNE with 4 solutions(!)
#-------------------------------------------------------------------------------

dimz <- cbind(x=c(1, 2), lam=c(3, 4), w=c(5, 6))
lagreq <- function(x, lam)
	c(	2*(x[1] - 2)*(x[2]-4)^4 + lam[1],
		2*(x[2] - 3)*x[1]^4 + lam[2]
		)
constr <- function(x)
	c(	x[1] + x[2] - 1,
		2*x[1] + x[2] - 2
		)		
jaclagreq <- function(x, lam)		
	rbind(	c(2*(x[2]-4)^4, 8*(x[1] - 2)*(x[2]-4)^3),
			c(8*(x[2] - 3)*x[1]^3, 2*x[1]^4)
			)
jacconstr <- function(x)
	rbind(	c(1, 1),
			c(2, 1)
			)
diaggrconstr <- function(x)
	diag(2)	

argH <- list(dimz=dimz, lagreq=lagreq, constr=constr)
argjac <- list(dimz=dimz, jaclagreq=jaclagreq, 
	jacconstr=jacconstr, diaggrconstr= diaggrconstr)
		
#base test
z <- rexp(6)
Htemplate(z, argH)
jacHtemplate(z, argjac)

potential.ce(1:10, 3, 3)
3*log(sum((1:10)^2)) - sum(log(4:10))
gradpotential.ce(1:10, 3, 3)
2*3/sum((1:10)^2)*1:10 - c(0, 0, 0, 1/4:10)

psi.ce(z, argH, 4)
gradpsi.ce(z, argH, argjac, 4)


#a simple test
x0 <- -c(2, 2)
z0 <- c(x0, 2, 2, pmax(10, 5-constr(x0) ) )

# jacHtemplate(z0, argjac)

GNE.ceq(z0, argH, argjac, silent=FALSE, control=list(trace=0))


x0 <- c(1, 0)
z0 <- c(x0, 2^9, 6, pmax(2, 2-constr(x0) ) )

# jacHtemplate(z0, argjac)

GNE.ceq(z0, argH, argjac, silent=FALSE, control=list(trace=1))




