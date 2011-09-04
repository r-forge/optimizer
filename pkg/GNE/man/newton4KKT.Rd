\name{NewtonKKT}
\alias{NewtonKKT}
\alias{newtonKKT}
\title{Newton type method}
\description{
		Newton type method for solving the extended KKT system of the GNEP.
	}
\usage{
NewtonKKT(xinit, method=c("Newton", "Levenberg-Marquardt","BFGS", "L-BFGS-B"), phi, jacphi, 
	control=list(), phi.arg=list(), jacphi.arg=list(), 
	extrapol=c("none", "RRE", "MPE", "SqRRE", "SqMPE"), bounds,
	silent=TRUE)
}
\arguments{
  \item{xinit}{initial point.}
  \item{method}{method name either \code{"Newton"}, \code{"Levenberg-Marquardt"}, \code{"BFGS"}, 
	\code{"L-BFGS-B"}. Default to \code{"Newton"}.}
  \item{phi}{the phi function.}
  \item{jacphi}{the Jacobian of the phi function.}
  \item{control}{a list of control parameters, see details.}
  \item{phi.arg}{further arguments to be passed to \code{phi}.}
  \item{jacphi.arg}{further arguments to be passed to \code{jacphi}.}
  \item{extrapol}{extrapolation name either \code{"none"}, \code{"Aitken"}, \code{"RRE"}, 
	\code{"MPE"}. Default to \code{"none"}.}
  \item{bounds}{a vector with bounds for L-BFGS-B method.}
  \item{silent}{a boolean to hide errors. Default to \code{TRUE}.}
}
\details{

The extended KKT system approach consists in solving the extended Karush-Kuhn-Tucker 
(KKT) system denoted by \eqn{\Phi(z)=0}, where eqn{z} is formed by the players strategy 
\eqn{x} and the Lagrange multiplier \eqn{\lambda}.
The root problem \eqn{\Phi(z)=0} is solved by an iterative scheme \eqn{z_{n+1} = z_n + d_n},
where the direction \eqn{d_n} is computed in three different ways:
\describe{
	\item{(a) Newton:}{The direction solves the system \eqn{Jac\Phi(x_n) d = - \Phi(x_n),} 
		generally called the Newton equation.}
	\item{(b) Levenberg-Marquardt:}{The direction solves the system 
		\deqn{ (Jac\Phi(x_n)^T Jac\Phi(x_n) + ||\Phi(x_n)|| I) d = - Jac\Phi(x_n)^T\Phi(x_n),}
		where \eqn{I} denotes the identity matrix.}
	\item{(c) quasi-Newton:}{quasi-Newton methods aim to solve an approximate version
		of the Newton equation \eqn{d = -\Phi(x_n) W_n} where \eqn{W_n} is computed
		by an iterative scheme. In the current implementation, \eqn{W_n} is updated
		by the BFGS method. The L-BFGS-B method is useful to add box constraints. }
}
Details on the methods can be found in Facchinei, Fischer & Piccialli (2009), \code{"Newton"}
corresponds to method 1 and \code{"Levenberg-Marquardt"} to method 3. We use the \code{R} implementation
of the quasi-Newton methods, see \code{\link{optim}}.


In addition to the Optimization approach, an extrapolation method may be used:
\describe{
	\item{- None}{no extrapolation.}
	\item{- RRE/MPE method}{Reduced Rank Extrapolation and Minimal Polynomial Extrapolation 
		method.}
	\item{- SqRRE/SqMPE method}{Squared version of RRE/MPE method.}	
}
See \code{\link{fixedpoint}} for details.


The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
	\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
	\item{\code{echo}}{A logical or an integer (0, 1, 2, 3) to print traces. 
		Default to \code{FALSE}, i.e. 0.}
}

}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{counts}}{A two-element integer vector giving the number of calls to \code{phi} and \code{jacphi} respectively.}			
		\item{\code{iter}}{The iteration number.}
		\item{\code{code}}{0 if convergence, 1 if \code{maxit} is reached, 10 
			if \code{tol} is not reached and 11 for both.}
	}
}
\references{

 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
 
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .


}
\seealso{

See \code{\link{fixedpoint}} and \code{\link{minGap}} for other approaches, and \code{\link{Phi}} and
\code{\link{JacPhi}} for template functions of \eqn{\Phi} and \eqn{Jac\Phi}.
}
\author{
 Christophe Dutang
}
\examples{

#-------------------------------------------------------------------------------
# (1) Example 5.1 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#Phi(z) function
phiex1 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]
	mu <- z[5]

	c(	2*(x-1) - s1 + mu,
		2*(y-1/2) - s2 + mu,
		min(x, s1),
		min(y, s2),
		min(1 - x - y, mu) )
}

#Jacobian of Phi(z)
jacphiex1 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]
	mu <- z[5]
	
	res <- matrix(0, 5, 5)
	res[1, ] <- c(2, 0, -1, 0, 1)
	res[2, ] <- c(0, 2, 0, -1, 1)
	res[3, ] <- c( 1*(x <= s1), 0, 1*(s1 <= x), 0, 0)
	res[4, ] <- c( 0, 1*(y <= s2), 0, 1*(s2 <= y), 0)
	res[5, ] <- c( -1*(1 - x - y <= mu), -1*(1 - x - y <= mu), 0, 0, -1*(mu <= 1 - x - y))
	res
}

#call, true value is (3/4, 1/4)
NewtonKKT(rep(0,5), "Newton", phiex1, jacphiex1, control=list(maxit=100, echo=TRUE))

NewtonKKT(rep(0,5), "Levenberg-Marquardt", phiex1, jacphiex1, control=list(maxit=100, echo=TRUE))

NewtonKKT(rep(-1,5), "Leven", phiex1, jacphiex1, control=list(maxit=100, echo=TRUE))


#-------------------------------------------------------------------------------
# (2) Example 5.2 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#constants
d <- 20
lambda <- 4
rho <- 1

#Phi(z) function  
phiex2 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]	
	c(	-(d - lambda - rho*(x+y) - rho*x) - s1,
		-(d - lambda - rho*(x+y) - rho*y) - s2,
		min(x, s1),
		min(y, s2)	)
}

#Jac Phi(z) function
jacphiex2 <- function(z)
{
	x <- z[1]
	y <- z[2]
	s1 <- z[3]
	s2 <- z[4]

	res <- matrix(0, 4, 4)
	res[1, ] <- c(2*rho, rho, -1, 0)
	res[2, ] <- c(rho, 2*rho, 0, -1)
	res[3, ] <- c( 1*(x <= s1), 0, 1*(s1 <= x), 0)
	res[4, ] <- c( 0, 1*(y <= s2), 0, 1*(s2 <= y))
	res
	
}

#call, true value is (16/3, 16/3) 

NewtonKKT(rep(0, 4), "Levenberg-Marquardt", phiex2, jacphiex2, control=list(maxit=200, echo=TRUE))

#-------------------------------------------------------------------------------
# (3) Example 5.3 of von Heusinger & Kanzow (2009)
#-------------------------------------------------------------------------------

#constants
cstC <- cbind(c(.1, .12, .15), c(.01, .05, .01))
cstU <- cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75))
cstK <- c(100, 100)
cstE <- c(.5, .25, .75)
cstD <- c(3, .01)

constr <- function(x, index)
	sum(cstU[, index] * cstE * x) - cstK[index]

gradL <- function(z, index)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	
	p1 <- -( cstD[1] - cstD[2]*sum(x) - cstC[index, 1] - 2*cstC[index, 2]*x[index] - cstD[2]*x[index])
	p2 <-  mu1*cstU[index,1]*cstE[index] + mu2*cstU[index,2]*cstE[index] 
	p1 + p2
}	

#phi function
phi <- function(z)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	
c(	gradL(z, 1),
	gradL(z, 2),
	gradL(z, 3),
	min(-constr(x, 1), mu1),
	min(-constr(x, 2), mu2) )
}	

#Jacobian of the phi function
jacphi <- function(z)
{
	x <- z[1:3]
	mu1 <- z[4]
	mu2 <- z[5]	
	A <- -constr(x, 1) <= mu1
	B <- -constr(x, 2) <= mu2
	
	res <- matrix(0, 5, 5)
	res[1, ] <- c(2*cstD[2]+2*cstC[1,2], cstD[2], cstD[2], cstU[1,1]*cstE[1], cstU[1,2]*cstE[1])
	res[2, ] <- c(cstD[2], 2*cstD[2]+2*cstC[2,2], cstD[2], cstU[2,1]*cstE[2], cstU[2,2]*cstE[2])
	res[3, ] <- c(cstD[2], cstD[2], 2*cstD[2]+2*cstC[3,2], cstU[3,1]*cstE[3], cstU[3,2]*cstE[3])
	res[4, ] <- c(-cstU[1,1]*cstE[1]*1*A, -cstU[2,1]*cstE[2]*1*A, -cstU[3,1]*cstE[3]*1*A, 1*(!A), 0)
	res[5, ] <- c(-cstU[1,2]*cstE[1]*1*B, -cstU[2,2]*cstE[2]*1*B, -cstU[3,2]*cstE[3]*1*B, 0, 1*(!B))

	res
}

#call, true value around (21.146, 16.027, 2.724)

NewtonKKT(rep(0,5), "Newton", phi, jacphi, control=list(maxit=2000, echo=FALSE))


NewtonKKT(rep(0,5), "Levenberg-Marquardt", phi, jacphi, control=list(maxit=2000, echo=FALSE))
	

}
\keyword{nonlinear}
\keyword{optimize}

