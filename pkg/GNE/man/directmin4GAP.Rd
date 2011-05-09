\name{algorithm for minimization}
\alias{minGap}
\alias{Gap}
\title{Algorithm for minimization of the reformulated problem}
\description{
		Algorithm for minimization of the reformulated problem using a gap function.
	}
\usage{
minGap(xinit, Gap, gradGap, method=c("BB","BFGS"), 
	control=list(), ...)


}
\arguments{
  \item{xinit}{initial point.}
  \item{Gap}{the gap/value function.}
  \item{gradGap}{the gradient of the gap/value function.}
  \item{method}{either \code{"BB"} or \code{BFGS}.}
  \item{control}{a list of control parameters, see details.}
  \item{\dots}{further arguments to be passed to \code{Gap} or \code{gradGap}.}
}
\details{

The gap function minimization consists in minimizing a gap function \eqn{min V(x)}. The function \code{minGap}
provides two optimization methods to solve this minimization problem.
\describe{
	\item{Barzilai-Borwein algorithm}{when \code{method = "BB"}, we use Barzilai-Borwein iterative scheme
	to find the minimum.}
	\item{Broyden-Fletcher-Goldfarb-Shanno algorithm}{when \code{method = "BFGS"}, we use the BFGS iterative
	scheme implemented in \code{R}, a quasi-Newton method with line search. }	
}
In the game theory literature, there are two main gap functions: the regularized 
Nikaido-Isoda (NI) function and the QVI regularized gap function. See von Heusinger & Kanzow (2009) 
for details on the NI function and Kubota & Fukushima (2009) for the QVI regularized gap function.


The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
	\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
	\item{\code{echo}}{A logical or an integer (0, 1, 2, 3) to print traces. 
		Default to \code{FALSE}, i.e. 0.}
	\item{\code{stepinit}}{Initial step size for the BB method (should be 
		small if gradient is ``big''). Default to 1.}
}


Note that the \code{Gap} function can return a numeric or a list with computation details. In the
latter case, the object return must be a list with the following components
\code{value}, \code{counts}, \code{iter}, see the example below.

}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{outer.counts}}{A two-element integer vector giving the number of 
			calls to \code{Gap} and \code{gradGap} respectively.}				
		\item{\code{outer.iter}}{The outer iteration number.}
		\item{\code{code}}{0 if convergence, 1 if \code{maxit} is reached, 10 
			if \code{tol} is not reached and 11 for both.}
		\item{\code{inner.iter.fn}, \code{inner.iter.gr}}{The iteration number
			for the function and the gradient when computing the gap function or
			its gradient (if appropriate).}	
		\item{\code{inner.counts.fn}, \code{inner.counts.gr}}{A two-element integer 
			vector giving the number of calls to the function and the gradient 
			when computing the gap function or its gradient (if appropriate).}	
	}
}
\references{
 A. von Heusinger (2009),
 \emph{Numerical Methods for the Solution of the Generalized Nash Equilibrium Problem},
 Ph. D. Thesis.
 
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .
 
 K. Kubota & M. Fukushima (2009), 
 \emph{Gap function approach to the generalized Nash Equilibrium problem},
 Journal of Optimization theory and applications.
}
\seealso{

See \code{\link{fixedpoint}} and \code{\link{NewtonKKT}} for other approaches.

}
\author{
 Christophe Dutang
}
\examples{

#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#Nikaido-Isoda function
NIF <- function(y,x, param)
	(x[1]-1)^2 - (y[1]-1)^2 + (x[2]-1/2)^2 - (y[2]-1/2)^2 - param/2*sum((x-y)^2)

#gradient of NIF
gradyNIF <- function(y, x, param)
c(	-2*(y[1]-1) + param*(x[1]-y[1]),
	-2*(y[2]-1/2) + param*(x[2]-y[2]) )	

#constraint function
constr <- function(y, x, param)
c(	-y[1],
	-y[2],
	y[1]+y[2]-1	)	

#gradient of constr
gradyconstr <- function(y, x, param)
rbind(
	c(-1, 0),
	c(0, -1),
	c(1, 1)
)

#starting value for the computation of yhat
xstart <- runif(2)/2

#yhat function
yhat <- function(x, alpha)
	constrOptim.nl(xstart, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))$par

#yhat function with full info
yhatfullres <- function(x, alpha)
	constrOptim.nl(xstart, fn=function(y, x, param) -NIF(y, x, param), 
		gr=function(y, x, param) -gradyNIF(y, x, param), 
		hin=function(y, x, param) -constr(y, x, param), 
		hin.jac=function(y, x, param) -gradyconstr(y, x, param), 
		x=x, param=alpha, control.outer=list(trace=FALSE))

#V hat function
V1alphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	p1 <- NIF(yalpha, x, alpha)	
	ybeta <- yhat(x, beta)
	p2 <- NIF(ybeta, x, beta)
	p1-p2
}

V1abfullres <- function(x, alpha, beta)
{
	yalpha <- yhatfullres(x, alpha)
	p1 <- NIF(yalpha$par, x, alpha)	
	ybeta <- yhatfullres(x, beta)
	p2 <- NIF(ybeta$par, x, beta)
	list(value=p1-p2, counts=yalpha$counts + ybeta$counts, iter=yalpha$outer.iterations + ybeta$outer.iterations)
}

#gradient of V1alphabeta
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

#call, true value is (3/4, 1/4)

#use BFGS method to minimize V hat function
minGap(c(0,0), V1alphabeta, GV1alphabeta, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))

#use BFGS method to minimize V hat function, with detailed info
minGap(c(0,0), V1abfullres, GV1abfullres, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))


#use BB method to minimize V hat function
minGap(c(0,0), V1abfullres, GV1abfullres, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))


#-------------------------------------------------------------------------------
# (2) Example 5.2 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

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


Valphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	ybeta <- yhat(x, beta)
	NIF(yalpha, x, alpha) - NIF(ybeta, x, beta)
}

GValphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	ybeta <- yhat(x, beta)
	gradxNIF(yalpha, x, alpha) - gradxNIF(ybeta, x, beta)
}


#call, true value is (16/3, 16/3) 

#use BB method to minimize V hat function via the Nikaido-Isoda function
minGap(c(1,1), Valphabeta, GValphabeta, method="BB", alpha=.02, beta=.05, control=list(tol=1e-8, echo=TRUE))



#test the QVI regularized gap function
yfunc <- function(x, alpha, F_x)
{
	if(missing(F_x))
		F_x <- sapply(1:2, function(i) gradtheta(x, i) )
	
	if(alpha != 0)
		z <- x - 1/alpha * F_x
	else
		stop("wrong alpha.")
	
	projector(z, g = function(y) constr(y), jacg = function(y) gradyconstr(y))
}


#RGP
reggapfunc <- function(x, alpha)
{
	F_x <- sapply(1:2, function(i) gradtheta(x, i) )
		
	y_x <- yfunc(x, alpha, F_x)
	as.numeric( crossprod(F_x, x - y_x) - alpha/2 * crossprod(x - y_x, x - y_x) )
}


#gradient of RGP w.r.t. x
gradRGP <- function(x, alpha)
{
	F_x <- sapply(1:2, function(i) gradtheta(x, i) )
	
	JacF_x <- t( sapply(1:2, function(j) sapply(1:2, function(i) grgrtheta(x, i, j) ) ) )
		
	y_x <- yfunc(x, alpha, F_x)

	as.numeric( -crossprod(JacF_x, y_x - x) + F_x - alpha*(x - y_x) )
}


grgrtheta <- function(z, i, j)
	rho*z[j] + rho*(i == j)

#with BB Algorithm
minGap(1:2, reggapfunc, gradRGP, method="BB", alpha=.1, control=list(tol=1e-8, echo=1, stepinit=1/100))

#with BFGS Algorithm
minGap(1:2, reggapfunc, gradRGP, method="BFGS", alpha=.1, control=list(tol=1e-8, echo=1))



#-------------------------------------------------------------------------------
# (3) Example 5.3 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

cstC <- cbind(c(.1, .12, .15), c(.01, .05, .01))
cstU <- cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75))
cstK <- c(100, 100)
cstE <- c(.5, .25, .75)
cstD <- c(3, .01)

#theta_j, index=j
theta <- function(z, index)
	- ( cstD[1] - cstD[2]*sum(z) - cstC[index, 1] - cstC[index, 2]*z[index] ) * z[index]

#grad theta_j w.r.t. z_j, index=j
gradtheta <- function(z, index)
	-( cstD[1] - cstD[2]*sum(z) - cstC[index, 1] - 2*cstC[index, 2]*z[index] - cstD[2]*z[index])

#grad theta_j w.r.t. z_-j, index=-j
grad2theta <- function(z, idj, id_j)
	cstD[2] * z[idj] 

#constraint functions
constr <- function(z, index)
	sum(cstU[, index] * cstE * z) - cstK[index]

gradconstr <- function(x, index)
	cstU[, index] * cstE * c(1,1,1)

fullconstr <- function(y, x, param)
	c(	constr(y, 1),
		constr(y, 2) )	

jacfullconstr <- function(y, x, param)
	rbind(	gradconstr(y, 1),
			gradconstr(y, 2) )	


#Nikaido-Isoda function
NIF <- function(y, x, param)
{
	p1 <- theta(x, 1) - theta(c(y[1], x[2], x[3]), 1) - param/2*(x[1]-y[1])^2
	p2 <- theta(x, 2) - theta(c(x[1], y[2], x[3]), 2) - param/2*(x[2]-y[2])^2
	p3 <- theta(x, 3) - theta(c(x[1], x[2], y[3]), 3) - param/2*(x[3]-y[3])^2
	p1 + p2 + p3	
}

#gradient w.r.t. y
gradyNIF <- function(y, x, param)
	c(- gradtheta(c(y[1], x[2], x[3]), 1) + param*(x[1]-y[1]),
	  - gradtheta(c(x[1], y[2], x[3]), 2) + param*(x[2]-y[2]),
	  - gradtheta(c(x[1], x[2], y[3]), 3) + param*(x[3]-y[3]) )
		
#gradient w.r.t. x		
gradxNIF <- function(y, x, param)
c(	
	gradtheta(x, 1) + grad2theta(x, 2, 1) - grad2theta(c(x[1], y[2], x[3]), 2, 1) + grad2theta(x, 3, 1) - grad2theta(c(x[1], x[2], y[3]), 3, 1) - param*(x[1] - y[1]),
	grad2theta(x, 1, 2) - grad2theta(c(y[1], x[2], x[3]), 1, 2) + gradtheta(x, 2) + grad2theta(x, 3, 2) - grad2theta(c(x[1], x[2], y[3]), 3, 2) - param*(x[2] - y[2]),
	grad2theta(x, 1, 3) - grad2theta(c(y[1], x[2], x[3]), 1, 3) + grad2theta(x, 2, 3) - grad2theta(c(x[1], y[2], x[3]), 2, 3) - gradtheta(x, 3) - param*(x[3] - y[3])
	)



startvalue <- 1:3


yhat <- function(x, alpha) 
{
	#random guess
	startvalue <- rejection(function(y) fullconstr(y,x,alpha), nvars=3, LB=0, UB=100)

	constrOptim.nl(par=startvalue, fn=function(y,z,param) -NIF(y,z,param),
		gr=function(y,z,param) -gradyNIF(y,z,param),
		hin=function(y,z,param) -fullconstr(y,z,param),
		hin.jac=function(y,z,param) -jacfullconstr(y,z,param),
		control.outer=list(trace=FALSE, eps=.1), 
		control.optim=list(trace=0, reltol=1e-6), z=x, param=alpha)$par
}

	
	
Valphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	p1 <- NIF(yalpha, x, alpha)	
	ybeta <- yhat(x, beta)
	p2 <- NIF(ybeta, x, beta)
	p1-p2
}
	
GValphabeta <- function(x, alpha, beta)
{
	yalpha <- yhat(x, alpha)
	p1 <- gradxNIF(yalpha, x, alpha)	
	ybeta <- yhat(x, beta)
	p2 <- gradxNIF(ybeta, x, beta)
	p1-p2
}

#call, true value around (21.146, 16.027, 2.724)
minGap(c(0,0,0), Valphabeta, GValphabeta, method="BB", alpha=.02, beta=.05, control=list(tol=1e-6, echo=TRUE))

minGap(c(0,0,0), Valphabeta, GValphabeta, method="BFGS", alpha=.02, beta=.05, control=list(tol=1e-6, echo=TRUE))


}
\keyword{nonlinear}
\keyword{optimize}

