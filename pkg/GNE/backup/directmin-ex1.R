library(GNE)

#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#

#constrained pb
NIF <- function(y,x, param)
{
#	cat("x", x, "y", y, "\n")
	(x[1]-1)^2 - (y[1]-1)^2 + (x[2]-1/2)^2 - (y[2]-1/2)^2 - param/2*sum((x-y)^2)
}

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
rbind(
	c(-1, 0),
	c(0, -1),
	c(1, 1)
)



start <- runif(1)
start <- c(start/3, start/3)
start <- c(runif(1)/2, runif(1)/2)


yhat <- function(x, alpha)
	constrOptim.nl(start, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))$par


yhatfullres <- function(x, alpha)
	constrOptim.nl(start, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))



yhatfullres(c(1/4,1/4), .01)

V1alphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	p1 <- NIF(x, yalpha, alpha)	
	ybeta <- yhat(x, beta)
	p2 <- NIF(x, ybeta, beta)
	p1-p2
}

V1abfullres <- function(x, alpha, beta)
{
	yalpha <- yhatfullres(x, alpha)
	p1 <- NIF(x, yalpha$par, alpha)	
	ybeta <- yhatfullres(x, beta)
	p2 <- NIF(x, ybeta$par, beta)
	list(value=p1-p2, counts=yalpha$counts + ybeta$counts, iter=yalpha$outer.iterations + ybeta$outer.iterations)
}

GV1alphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	g1 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - alpha*(x-yalpha)
	ybeta <- yhat(x, beta)
	g2 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - beta*(x-ybeta)
	g1-g2
}
GV1abfullres <- function(x, alpha, beta)
{
	yalpha <- yhatfullres(x, alpha)
	g1 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - alpha*(x-yalpha$par)
	ybeta <- yhatfullres(x, beta)
	g2 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - beta*(x-ybeta$par)
	list(value=g1-g2, counts=yalpha$counts + ybeta$counts, iter=yalpha$outer.iterations + ybeta$outer.iterations)
}


V1alphabeta(c(1/4,3/4), .02, .05)
V1abfullres(c(1/4,3/4), .02, .05)

GV1alphabeta(c(1/4,3/4), .02, .05)
GV1abfullres(c(1/4,3/4), .02, .05)

minGap(c(0,0), V1alphabeta, GV1alphabeta, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))

minGap(c(0,0), V1abfullres, GV1abfullres, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))



minGap(c(0,0), V1abfullres, GV1abfullres, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))


