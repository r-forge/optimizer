
require(NE)	

#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#

#as function of y for a given x and param=alpha

NIF <- function(y, x, param) 
	(x[1]-1)^2 - (y[1]-1)^2 + (x[2]-1/2)^2 - (y[2]-1/2)^2 - param/2*sum((x-y)^2)

#NIF(runif(2), runif(2), .1)

gradyNIF <- function(y, x, param)
{
c(	-2*(y[1]-1)-param*(x[1]-y[1]),
	-2*(y[2]-1/2)-param*(x[2]-y[2]) )	
}

constr <- function(y, x, param)
c(	-y[1],
	-y[2],
	y[1]+y[2]-1	)	


gradyconstr <- function(y, x, param)
cbind(
	c(-1, 0, 1),
	c(0, -1, 1)		
)


yhat <- function(x, alpha)
	constrOptim.nl(runif(2)/2, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))$par


		
yhatfullres <- function(x, alpha)
{
	res <- constrOptim.nl(runif(2)/2, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE),
		control.optim=list(reltol=1e-6, trace=1))
		
	list(value=res$value, counts=res$counts, iter=res$outer.iterations, par=res$par)
}		
		

Vfunc <- function(x, alpha)
	NIF(yhat(x, alpha), x, alpha)

Vfuncfullres <- function(x, alpha)
{
	yalpha <- yhatfullres(x, alpha)
	
	list(value=NIF(yalpha$par, x, alpha), counts=yalpha$counts, iter=yalpha$outer.iterations)
}



fixedpointrelaxation(c(0,0), method="UR", yhat, purestep, Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpointrelaxation(c(0,0), method="UR", yhatfullres, purestep, Vfuncfullres, control=list(echo=TRUE), alpha=.1)



fixedpointrelaxation(c(0,0), method="vH", yhat, Vfunc=Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpointrelaxation(c(0,0), method="vH", yhatfullres, Vfunc= Vfuncfullres, control=list(echo=TRUE), alpha=.1)
