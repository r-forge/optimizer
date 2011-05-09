\name{fixedpoint}
\alias{fixedpoint}
\title{Fixed point method to solve GNE problem.}
\description{
		Fixed point method to solve GNE problem.
	}
\usage{
fixedpoint(xinit, method=c("pure","UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"), yfunc, stepfunc, Vfunc, 
	control=list(), ...)

}
\arguments{
  \item{xinit}{initial point.}
  \item{method}{either \code{"pure"}, \code{"UR"}, \code{"vH"}, \code{"RRE"}, \code{"MPE"}, 
	\code{"SqRRE"} or \code{"SqMPE"} method, see details.}
  \item{yfunc}{the y function.}
  \item{stepfunc}{the step function, only needed when \code{method="UR"}.}
  \item{Vfunc}{the value function used to compute the merit function.}
  \item{control}{a list of control parameters, see details.}
  \item{\dots}{further arguments to be passed to \code{yfunc}, \code{stepfunc} or \code{Vfunc}.}
}
\details{

The fixed point approach consists in solving equation \eqn{y(x)=x}.

\describe{
	\item{(a) Crude or pure fixed point method:}{
		It simply consists in iterations \eqn{x_{n+1} = y(x_n)}.}
	\item{(b) Polynomial methods:}{
			\describe{
				\item{- relaxation algorithm (linear extrapolation):}{
					The next iterate is computed as \deqn{x_{n+1} = (1-\alpha_n) x_n + \alpha_n y(x_n).}
					The step \eqn{\alpha_n} can be computed in different ways: constant, decreasing
					serie or a line search method. In the literature of game theory, the decreasing serie
					refers to the method of Ursayev and Rubinstein (\code{method="UR"}) while the line search
					method refers to the method of von Heusinger (\code{method="vH"}). Note that the constant
					step can be done using the UR method.}
				\item{- RRE and MPE method:}{
					Reduced Rank Extrapolation and Minimal Polynomial Extrapolation 
					methods are polynomial extrapolation methods, where the monomials are functional 
					``powers'' of the y function, i.e. function composition of y. Of order 1, RRE and MPE
					consists of \deqn{x_{n+1} = x_n + t_n (y(x_n) - x_n),} 
					where \eqn{t_n} equals to
					\eqn{<v_n, r_n> / <v_n, v_n>} for RRE1 and \eqn{<r_n, r_n> / <v_n, r_n>} for MPE1, where
					\eqn{r_n =y(x_n) - x_n } and \eqn{v_n = y(y(x_n)) - 2y(x_n) + x_n}. 
					To use RRE/MPE methods, set \code{method = "RRE"} or \code{method = "MPE"}.}	
				\item{- squaring method:}{
					It consists in using an extrapolation method (such as RRE and MPE)
					after two iteration of the linear extrapolation, i.e. 
					\deqn{x_{n+1} = x_n -2 t_n r_n + t_n^2 v_n.} The squared version of RRE/MPE methods are
					available via setting \code{method = "SqRRE"} or \code{method = "SqMPE"}.}	
			}
		}
	\item{(c) Epsilon algorithms:}{Not implemented.}
}
For details on fixed point methods, see Varadhan & Roland (2004).


The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
	\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
	\item{\code{echo}}{A logical or an integer (0, 1, 2, 3) to print traces. 
		Default to \code{FALSE}, i.e. 0.}
	\item{\code{sigma, beta}}{parameters for von Heusinger algorithm. 
		Default to 9/10 and 1/2 respectively.}
}

Note that the \code{yfunc}, \code{Vfunc} functions can return a numeric or a list with computation details. In the
latter case, the object return must be a list with the following components
\code{value}, \code{counts}, \code{iter}, see the example below.


}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{counts}}{A two-element integer vector giving the number of calls to \code{yfunc} and \code{Vfunc}/\code{stepfunc} respectively.}			
		\item{\code{iter}}{The iteration number.}
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
 
 S. Uryasev & R.Y. Rubinstein (1994), 
 \emph{On relaxation algorithms in computation of noncooperative equilibria}, 
 IEEE Transactions on Automatic Control.
 
 R. Varadhan & C. Roland (2004),
 \emph{Squared Extrapolation Methods (SQUAREM): A New Class of Simple and Efficient Numerical 
 Schemes for Accelerating the Convergence of the EM Algorithm},
 Johns Hopkins University, Dept. of Biostatistics Working Papers.
 
}
\seealso{

See \code{\link{NewtonKKT}} and \code{\link{minGap}} for other approaches, and \code{\link{stepfunc}} for step functions.

}
\author{
 Christophe Dutang
}
\examples{

#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#as function of y for a given x and param=alpha

NIF <- function(y, x, param) 
	(x[1]-1)^2 - (y[1]-1)^2 + (x[2]-1/2)^2 - (y[2]-1/2)^2 - param/2*sum((x-y)^2)

#NIF(runif(2), runif(2), .1)

gradyNIF <- function(y, x, param)
{
c(	-2*(y[1]-1) + param*(x[1]-y[1]),
	-2*(y[2]-1/2) + param*(x[2]-y[2]) )	
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


#call, true value is (3/4, 1/4)

fixedpoint(c(0,0), method="pure", yhat, Vfunc=Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpoint(c(0,0), method="pure", yhatfullres, Vfunc= Vfuncfullres, control=list(echo=TRUE), alpha=.1)


fixedpoint(c(0,0), method="UR", yhat, purestep, Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpoint(c(0,0), method="UR", yhatfullres, purestep, Vfuncfullres, control=list(echo=TRUE), alpha=.1)



fixedpoint(c(0,0), method="vH", yhat, Vfunc=Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpoint(c(0,0), method="vH", yhatfullres, Vfunc= Vfuncfullres, control=list(echo=TRUE), alpha=.1)


fixedpoint(c(0,0), method="SqRRE", yhat, Vfunc=Vfunc, control=list(echo=TRUE), alpha=.1)

fixedpoint(c(0,0), method="SqRRE", yhatfullres, Vfunc= Vfuncfullres, control=list(echo=TRUE), alpha=.1)
	

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



Vfunc <- function(x, alpha)
	NIF(yhat(x, alpha), x, alpha)

#call, true value is (16/3, 16/3) 

fixedpoint(c(0,0), method="pure", yhat, decrstep20, Vfunc, control=list(echo=FALSE), alpha=.01)
fixedpoint(c(0,0), method="pure", yhat, decrstep20, Vfunc, control=list(echo=TRUE), alpha=.01)
fixedpoint(c(0,0), method="pure", yhat, decrstep20, Vfunc, control=list(echo=2), alpha=.01)
fixedpoint(c(0,0), method="pure", yhat, decrstep20, Vfunc, control=list(echo=3), alpha=.01)


fixedpoint(c(0,0), method="UR", yhat, decrstep20, Vfunc, control=list(echo=TRUE), alpha=.01)


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
	res <- c(	- gradtheta(c(y[1], x[2], x[3]), 1) + param*(x[1]-y[1]),
		- gradtheta(c(x[1], y[2], x[3]), 2) + param*(x[2]-y[2]),
		- gradtheta(c(x[1], x[2], y[3]), 3) + param*(x[3]-y[3]) )
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

	
	
Valpha <- function(x, alpha, mytrace)
{
	yalpha <- yhat(x, alpha, mytrace)
	NIF(x, yalpha, alpha)	
}


#call, true value around (21.146, 16.027, 2.724)
fixedpoint(c(0, 0, 0), method="UR", yhat, decrstep10, Valpha, control=list(echo=TRUE), alpha=.02, mytrace=FALSE)




}
\keyword{nonlinear}
\keyword{optimize}

