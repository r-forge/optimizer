
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


#first call
xk <- c(3.7826847, 0.5112527)
zk <- c(xk, 10, 10, pmax(10, 5-constr(xk) ) )
Htemplate(zk, argH)
psi.ce(zk, argH, length(zk)/2)

GNE.ceq(zk, argH, argjac, silent=FALSE, control=list(trace=3, btol=1e-2))



zeta <- 4	
zk <- c(3.7826847, 0.5112527, 0.9206536, 0.2450269, 1.0918416, 1.7969769)	
dk <- ceq.direction(zk, argH, argjac, 1/2)
slopek <- as.numeric( crossprod( gradpsi(zk, argH, argjac, zeta) , dk ) )
control <- list(ftol=1e-6, xtol=1e-5, btol=1e-2, maxit=100, trace=0, sigma=1/2, 
				echofile=NULL, delta=1)

psi.ce(zk+dk, argH=argH, zeta=zeta)


linesearch.geom(zk, dk, slopek, control, psi, argH=argH, zeta=zeta)	

linesearch.geom(zk, dk, slopek, control, function(z) as.numeric(crossprod(z)/2)	)



