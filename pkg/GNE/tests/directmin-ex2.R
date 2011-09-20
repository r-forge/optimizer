library(GNE)

#constants
d <- 20
lambda <- 4
rho <- 1
start <- runif(2)

#theta_i function
theta <- function(z, index)
	- (d-lambda-rho*sum(z))*z[index]
	
#grad theta_i w.r.t. to z_i
gradtheta <- function(z, index)
	-( d-lambda-rho*sum(z)-rho*z[index] )
	
#grad theta_i w.r.t. to z_-i
grad2theta <- function(z, index)
	-rho*z[index]
	
#Nikaido-Isoda function
NIF <- function(y, x, param)
	theta(x, 1) - theta( c(y[1], x[2]), 1) + theta(x, 2) - theta( c(x[1], y[2]), 2) - param/2*sum((x-y)^2)

#gradient w.r.t. y
gradyNIF <- function(y, x, param)
c(	-gradtheta( c(y[1], x[2]), 1) + param*(x[1]-y[1]),
	-gradtheta( c(x[1], y[2]), 2) + param*(x[2]-y[2]) )

#gradient w.r.t. x
gradxNIF <- function(y, x, param)
c(	gradtheta(x, 1) + grad2theta(x, 2) - grad2theta(c(x[1], y[2]), 2) - param*(x[1]-y[1]),
	grad2theta(x, 1) - grad2theta(c(y[1], x[2]), 1) + gradtheta(x, 2) - param*(x[2]-y[2]) )

#constraint function
constr <- function(y, x, param)
c(	-y[1],
	-y[2]	)

#gradient of constr
gradyconstr <- function(y, x, param)
cbind(
	c(-1, 0),
	c(0, -1)		
)


yhat <- function(x, alpha)
{
	xstart <- x+rexp(2)
	res <- constrOptim.nl(xstart, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))
	res$par
}

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


Valphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	ybeta <- yhat(x, beta)
	NIF(yalpha, x, alpha) - NIF(ybeta, x, beta)
}

Vabfullres <- function(x, alpha, beta)
{
	yalpha <- yhatfullres(x, alpha)
	p1 <- NIF(yalpha$par, x, alpha)	
	ybeta <- yhatfullres(x, beta)
	p2 <- NIF(ybeta$par, x, beta)
	list(value=p1-p2, counts=yalpha$counts + ybeta$counts, iter=yalpha$outer.iterations + ybeta$outer.iterations)
}

GValphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	ybeta <- yhat(x, beta)
	gradxNIF(yalpha, x, alpha) - gradxNIF(ybeta, x, beta)
}


GVabfullres <- function(x, alpha, beta)
{
	yalpha <- yhatfullres(x, alpha)
	g1 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - alpha*(x-yalpha$par)
	ybeta <- yhatfullres(x, beta)
	g2 <- c( 2*(x[1]-1), 2*(x[2]-1/2) ) - beta*(x-ybeta$par)
	list(value=g1-g2, counts=yalpha$counts + ybeta$counts, iter=yalpha$outer.iterations + ybeta$outer.iterations)
}


#call, true value is (16/3, 16/3) 

#use BB method to minimize V hat function via the Nikaido-Isoda function
minGap(c(10,20), Valphabeta, GValphabeta, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))

# with detailed info
minGap(c(10,20), Vabfullres, GVabfullres, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=FALSE))

gapBB <- minGap(c(0,0), Vabfullres, GVabfullres, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=FALSE))

gapBB$inner.counts.fn + gapBB$inner.counts.gr

gapBFGS <- minGap(c(0,0), Vabfullres, GVabfullres, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-8, echo=FALSE))

gapBFGS$inner.counts.fn + gapBFGS$inner.counts.gr


#test the QVI regularized gap function
yfunc <- function(x, alpha, grad_x)
{


	if(missing(grad_x))
		grad_x <- sapply(1:2, function(i) gradtheta(x, i) )

	
	if(alpha != 0)
		z <- x - 1/alpha * grad_x
	else
		stop("wrong alpha.")
	
	projector(z, 
		g = function(y, ...) constr(y, ...),
		jacg = function(y, ...) gradyconstr(y, ...))
}


#RGP
reggapfunc <- function(x, alpha)
{
	grad_x <- sapply(1:2, function(i) gradtheta(x, i) )
		
	y_x <- yfunc(x, alpha, grad_x)
	as.numeric( crossprod(grad_x, x - y_x) - alpha/2 * crossprod(x - y_x, x - y_x) )
}

#gradient of RGP w.r.t. x
gradRGP <- function(x, alpha)
{
	F_x <- sapply(1:2, function(i) gradtheta(x, i) )
	
	JacF_x <- sapply(1:2, function(j) 
		sapply(1:2, function(i) grgrtheta(x, i, j) )
		)
	
	y_x <- yfunc(x, alpha, F_x)
	

	as.numeric( -JacF_x %*% (y_x - x) + F_x - alpha*(x - y_x) )
}


grgrtheta <- function(z, i, j)
	rho*z[j] + rho*(i == j)

yfunc(1:2, .1)
yfunc(1:2, .1)


reggapfunc(c(1,2), .1)
reggapfunc(c(1,2), 1)
gradRGP(c(1,2), .1)

c(grgrtheta(1:2, 1, 1), grgrtheta(1:2, 1, 2))

minGap(1:2, reggapfunc, gradRGP, method="BB", alpha=.1, control=list(tol=1e-8, echo=1, stepinit=1/100))

minGap(1:2, reggapfunc, gradRGP, method="BFGS", alpha=.1, control=list(tol=1e-8, echo=1))

