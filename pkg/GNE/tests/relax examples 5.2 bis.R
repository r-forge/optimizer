require(GNE)

#-------------------------------------------------------------------------------
# (2) Example 5.2 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

d <- 20
lambda <- 4
rho <- 1
start <- runif(2)

theta <- function(z, index)
	- (d-lambda-rho*sum(z))*z[index]
	

gradtheta <- function(z, index)
	-( d-lambda-rho*sum(z)-rho*z[index] )
	

NIF <- function(y, x, param)
	theta(x, 1) - theta( c(y[1], x[2]), 1) + theta(x, 2) - theta( c(x[1], y[2]), 2) - param/2*sum((x-y)^2)

gradyNIF <- function(y, x, param)
c(	-gradtheta( c(y[1], x[2]), 1) + param*(x[1]-y[1]),
	-gradtheta( c(x[1], y[2]), 2) + param*(x[2]-y[2]) )



constr <- function(y, x, param)
c(	-y[1],
	-y[2]	)


gradyconstr <- function(y, x, param)
cbind(
	c(-1, 0),
	c(0, -1)		
)




yhat <- function(x, alpha)
	constrOptim.nl(x+rexp(2), fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))$par


yhatfullres <- function(x, alpha)
{
	xstart <- x+rexp(2)
	res <- constrOptim.nl(xstart, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))
	res
}



Vfunc <- function(x, alpha)
	NIF(yhat(x, alpha), x, alpha)
	
Vfullres <- function(x, alpha)
{
	yalpha <- yhatfullres(x, alpha)
	list(value=NIF(yalpha$par, x, alpha), counts=yalpha$counts, iter=yalpha$outer.iterations)
}	
	

#call, true value is (16/3, 16/3) 

fixedpoint(c(0,0), method="pure", yhat, decrstep20, Vfunc, control=list(echo=FALSE), alpha=.01)

purFP <- fixedpoint(c(0,0), method="pure", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

purFP$inner.counts.fn + purFP$inner.counts.vh


fixedpoint(c(0,0), method="UR", yhat, decrstep20, Vfunc, control=list(echo=TRUE), alpha=.01)

decr20 <- fixedpoint(c(0,0), method="UR", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

decr20$inner.counts.fn + decr20$inner.counts.vh


fixedpoint(c(0,0), method="UR", yhat, decrstep20, Vfunc, control=list(echo=2), alpha=.01)

fixedpoint(c(0,0), method="UR", yhat, decrstep20, Vfunc, control=list(echo=3), alpha=.01)


fixedpoint(c(0,0), method="vH", yhat, decrstep20, Vfunc, control=list(echo=TRUE, sigma=1/2), alpha=.01)

linesearch <- fixedpoint(c(0,0), method="vH", yhatfullres, decrstep20, Vfullres, control=list(echo=TRUE, sigma=1/2), alpha=.01)

linesearch$inner.counts.fn + linesearch$inner.counts.vh


fixedpoint(c(10,10), method="vH", yhat, decrstep20, Vfunc, control=list(echo=TRUE, sigma=9/10, beta=2/3), alpha=.01)




RREFP <- fixedpoint(c(0,0), method="RRE", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

RREFP$inner.counts.fn + RREFP$inner.counts.vh

MPEFP <- fixedpoint(c(0,0), method="MPE", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

MPEFP$inner.counts.fn + MPEFP$inner.counts.vh


SqRREFP <- fixedpoint(c(0,0), method="SqRRE", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

SqRREFP$inner.counts.fn + SqRREFP$inner.counts.vh



SqMPEFP <- fixedpoint(c(0,0), method="SqMPE", yhatfullres, decrstep20, Vfullres, control=list(echo=FALSE), alpha=.01)

SqMPEFP$inner.counts.fn + SqMPEFP$inner.counts.vh

