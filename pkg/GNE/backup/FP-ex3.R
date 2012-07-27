library(GNE)

# (3) Example 5.3 of von Heusinger & Kanzow (2009)
#

cstC <- cbind(c(.1, .12, .15), c(.01, .05, .01))
cstU <- cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75))
cstK <- c(100, 100)
cstE <- c(.5, .25, .75)
cstD <- c(3, .01)

#theta_j, index=j
theta <- function(z, index)
{
	- ( cstD[1] - cstD[2]*sum(z) - cstC[index, 1] - cstC[index, 2]*z[index] ) * z[index]
}

#grad theta_j w.r.t. z_j, index=j
gradtheta <- function(z, index)
	-( cstD[1] - cstD[2]*sum(z) - cstC[index, 1] - 2*cstC[index, 2]*z[index] - cstD[2]*z[index])

#grad theta_j w.r.t. z_-j, index=-j
grad2theta <- function(z, idj, id_j)
	cstD[2] * z[idj] 


constr <- function(z, index)
	sum(cstU[, index] * cstE * z) - cstK[index]


gradconstr <- function(x, index)
	cstU[, index] * cstE * c(1,1,1)

#constraint function for y : g(y)
fullconstr <- function(y, x, param)
	c(	constr(y, 1),
		constr(y, 2) )	
		
fullconstr(c(21.1, 16.0, 2.7))
		
#jacobian constraint function for y : Grad g(y)
jacfullconstr <- function(y, x, param)
	rbind(	gradconstr(y, 1),
			gradconstr(y, 2) )	

jacfullconstr(c(21.1, 16.0, 2.7))


NIF <- function(y, x, param)
{

	p1 <- theta(x, 1) - theta(c(y[1], x[2], x[3]), 1) - param/2*(x[1]-y[1])^2
	p2 <- theta(x, 2) - theta(c(x[1], y[2], x[3]), 2) - param/2*(x[2]-y[2])^2
	p3 <- theta(x, 3) - theta(c(x[1], x[2], y[3]), 3) - param/2*(x[3]-y[3])^2
	p1 + p2 + p3	
}


gradyNIF <- function(y, x, param)
{
	res <- c(	- gradtheta(c(y[1], x[2], x[3]), 1) - param*(x[1]-y[1]),
		- gradtheta(c(x[1], y[2], x[3]), 2) - param*(x[2]-y[2]),
		- gradtheta(c(x[1], x[2], y[3]), 3) - param*(x[3]-y[3]) )
	res	
}	



yhat <- function(x, alpha, mytrace, tol=1e-6) 
{
	constrOptim.nl(par=x, fn=function(y,z,param) -NIF(y,z,param),
		gr=function(y,z,param) -gradyNIF(y,z,param),
		hin=function(y,z,param) -fullconstr(y,z,param),
		hin.jac=function(y,z,param) -jacfullconstr(y,z,param),
		control.outer=list(trace=FALSE, eps=.1), 
		control.optim=list(trace=0, reltol=1e-6), z=x, param=alpha)$par
}

yhat(rexp(3), .02)		
	
	
Valpha <- function(x, alpha, mytrace)
{
	yalpha <- yhat(x, alpha, mytrace)
	NIF(x, yalpha, alpha)	
}


#fixedpoint(c(0, 0, 0), method="UR", yhat, decrstep10, control=list(echo=TRUE), alpha=.02, mytrace=FALSE)



#fixedpointrelaxation(c(0, 0, 0), method="vH", yhat, Vfunc=Valpha, control=list(echo=TRUE, beta=1/20), alpha=.01, mytrace=FALSE)

