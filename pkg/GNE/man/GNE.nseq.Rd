\name{GNE.nseq}
\alias{GNE.nseq}
\title{Non smooth equation reformulation of the GNE problem.}
\description{
		Non smooth equation reformulation via the extended KKT system of the GNE problem.
	}
\usage{
GNE.nseq(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), method="Newton", control=list(), ...)
bench.GNE.nseq(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), echo=FALSE, ...)
}
\arguments{
  \item{xinit}{initial point.}
  \item{Phi}{the phi function.}
  \item{jacPhi}{the Jacobian of the phi function.}
  \item{argPhi}{further arguments to be passed to \code{phi}.}
  \item{argjac}{further arguments to be passed to \code{jacphi}.}
  \item{method}{a character string specifying the method \code{"Newton"},
	\code{"Broyden"}, \code{"Levenberg-Marquardt"}. }
  \item{control}{a list with control parameters.}
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions \code{phi} and \code{jacphi}.}
  \item{echo}{a logical to get some traces. Default to \code{FALSE}.}	
}
\details{
\code{GNE.nseq} solves the GNE problem via a non smooth reformulation of the KKT system. 
\code{bench.GNE.nseq} carries out a benchmark of the computation methods for a given initial point.

This approach consists in solving the extended Karush-Kuhn-Tucker 
(KKT) system denoted by \eqn{\Phi(z)=0}, where eqn{z} is formed by the players strategy 
\eqn{x} and the Lagrange multiplier \eqn{\lambda}.
The root problem \eqn{\Phi(z)=0} is solved by an iterative scheme \eqn{z_{n+1} = z_n + d_n},
where the direction \eqn{d_n} is computed in three different ways:
\describe{
	\item{(a) Newton:}{The direction solves the system \eqn{Jac\Phi(x_n) d = - \Phi(x_n)}{Jac Phi(x_n) d = - Phi(x_n)}, 
		generally called the Newton equation.}
	\item{(b) Broyden:}{It is a quasi-Newton method aiming to solve an approximate version
		of the Newton equation \eqn{d = -\Phi(x_n) W_n}{d = -Phi(x_n) W_n} where \eqn{W_n} is computed
		by an iterative scheme. In the current implementation, \eqn{W_n} is updated
		by the Broyden method. }
	\item{(c) Levenberg-Marquardt:}{The direction solves the system 
		\deqn{ (Jac\Phi(x_n)^T Jac\Phi(x_n) + ||\Phi(x_n)||^\delta I) d = - Jac\Phi(x_n)^T\Phi(x_n)
			}{ (Jac Phi(x_n)^T Jac Phi(x_n) + ||Phi(x_n)||^delta I) d = - Jac Phi(x_n)^T Phi(x_n),}
		where \eqn{I} denotes the identity matrix.}
}
Details on the methods can be found in Facchinei, Fischer & Piccialli (2009), \code{"Newton"}
corresponds to method 1 and \code{"Levenberg-Marquardt"} to method 3. For the Broyden, see 
Dennis & Moree (1977).

In addition to the computation method, a globalization scheme can be choosed using the \code{global}
argument. Available schemes are 
\describe{
	\item{(1) Line search:}{ if \code{global} is set to \code{"qline"} or \code{"gline"}, a line search
		is used with the merit function being half of the L2 norm of \eqn{Phi}, respectively with a
		quadratic or a geometric implementation.}
	\item{(2) Trust region:}{ if \code{global} is set to \code{"dbldog"} or \code{"pwldog"}, a trust
		region is used respectively with a double dogleg or a Powell (simple) dogleg implementation.}
	\item{(3) None:}{ if  \code{global} is set to \code{"none"}, no globalization is done. }
}	
The default value of \code{global} is \code{"gline"}. 

The implementation relies heavily on the 
\code{nleqslv} function of the package of the same name. So full details on the control parameters are
to be found in the help page of this function. We briefly recall here the main parameters.
The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{xtol}}{The relative steplength tolerance.
	When the relative  steplength of all scaled x values is smaller than this value
	convergence is declared. The default value is \eqn{10^{-8}}{1e-8}.
	}
	\item{\code{ftol}}{The function value tolerance.
	Convergence is declared when the largest absolute function value is smaller than \code{ftol}.
	The	default value is \eqn{10^{-8}}{1e-8}.
	}
	\item{\code{maxit}}{The maximum number of major iterations. The default value is 150 if a 
	global strategy has been specified.}
	\item{\code{trace}}{Non-negative integer. A value of 1 will give a detailed report of the
	progress of the iteration.}
}

}
\value{
\code{GNE.nseq} returns a list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{value}}{The value of the merit function.}
		\item{\code{counts}}{A two-element integer vector giving the number of calls to 
			\code{phi} and \code{jacphi} respectively.}			
		\item{\code{iter}}{The outer iteration number.}
		\item{\code{code}}{
			         The values returned are
         \describe{
			\item{\code{1}}{Function criterion is near zero.
			Convergence of function values has been achieved.}
			\item{\code{2}}{x-values within tolerance. This means that the relative distance between two
			consecutive x-values is smaller than \code{xtol}.}
			\item{\code{3}}{No better point found.
			This means that the algorithm has stalled and cannot find an acceptable new point.
			This may or may not indicate acceptably small function values.}
			\item{\code{4}}{Iteration limit \code{maxit} exceeded.}
			\item{\code{5}}{Jacobian is too ill-conditioned.}
			\item{\code{6}}{Jacobian is singular.}
			\item{\code{100}}{an error in the execution.}
			}
		}
		\item{\code{message}}{a string describing the termination code}	
	}
\code{bench.GNE.nseq} returns a list with components:
	\describe{
		\item{\code{compres}}{a data.frame summarizing the different computations.}
		\item{\code{reslist}}{a list with the different results from \code{GNE.nseq}.}
	} 	

}
\references{

 J. E. Dennis and J. J. Moree (1977), 
 \emph{Quasi-Newton methods, Motivation and Theory},
 SIAM review.
 
 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
 
 B. Hasselman (2011), 
 \emph{nleqslv: Solve systems of non linear equations},
 R package.

 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem 
 using Nikaido-Isoda-type functions},
 Comput Optim Appl .


}
\seealso{

See \code{\link{GNE.fpeq}} and \code{\link{GNE.min}} 
for other approaches; \code{\link{Phi}} and
\code{\link{JacPhi}} for template functions of \eqn{\Phi} and \eqn{Jac\Phi} and
\code{\link{complementarity}} for complementarity functions.

See \code{\link{nleqslv}} for the optimization routine.
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

#true value is (3/4, 1/4, 0, 0, 1/2)

z0 <- rep(0, 5)
GNE.nseq(z0, phiex1, jacphiex1, method="Newton")

GNE.nseq(z0, phiex1, jacphiex1, method="Broyden")


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

#call, true value is (16/3, 16/3, 0, 0) 

z0 <- rep(0, 4)
GNE.nseq(z0, phiex2, jacphiex2, method="Newton")

GNE.nseq(z0, phiex2, jacphiex2, method="Broyden")


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
phiex3 <- function(z)
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
jacphiex3 <- function(z)
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

#call, true value around (21.146, 16.027, 2.724, 0.574, 0.000)

z0 <- rep(0, 5)
GNE.nseq(z0, phiex3, jacphiex3, method="Newton")

GNE.nseq(z0, phiex3, jacphiex3, method="Broyden")


#-------------------------------------------------------------------------------
# (4) Example of GNE with 4 solutions(!)
#-------------------------------------------------------------------------------



F <- function(z, phi=phiMin, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	# cat("x", x, "\n")
	c(	2*(x[1] - 2)*(x[2]-4)^4 + lambda[1],
		2*(x[2] - 3)*x[1]^4 + lambda[2],
		phi(lambda[1], 1-sum(x), ...),
		phi(lambda[2], 2-2*x[1]-x[2], ...)
	)
}


JacF <- function(z, gphia, gphib, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	idga1 <- gphia(lambda[1], 1- sum(x), ...)
	idgb1 <- gphib(lambda[1], 1- sum(x), ...)
	idga2 <- gphia(lambda[2], 2-2*x[1]-x[2], ...)
	idgb2 <- gphib(lambda[2], 2-2*x[1]-x[2], ...)

rbind(
c(2*(x[2]-4)^4, 8*(x[1] - 2)*(x[2]-4)^3, 1, 0),
c(8*(x[2] - 3)*x[1]^3, 2*x[1]^4, 0, 1),
c(-idgb1, -idgb1, idga1, 0),
c(-2*idgb2, -idgb2, 0, idga2)
)	
}


#list of true GNEs
trueGNE <- rbind(c(2, -2, 0, 5*2^5),
	c(-2, 3, 8, 0),
	c(0, 1, 4*3^4, 0),
	c(1, 0, 2^9, 6))
colnames(trueGNE) <- c("x1", "x2", "lam1", "lam2")
rownames(trueGNE) <- 1:4

z0 <- c(10, 10, 1, 1)

GNE.nseq(z0, F, JacF, list(phi=phiMin), list(gphia= GrAphiMin, gphib= GrBphiMin), method="Newton")

GNE.nseq(z0, F, JacF, list(phi= phiFB), list(gphia= GrAphiFB, gphib= GrBphiFB), method="Newton")


#random initial points
n <- 20
set.seed(1234)
initpt <- cbind(runif(n, -10, 10), runif(n, -10, 10), 1, 1)

NewLnsrch <- function(i, echo=FALSE)
{ 
	if(echo)
		cat("______", initpt[i, ], "\n")
	res <- GNE.nseq(initpt[i, ], F, JacF, list(phi=phiMin), list(gphia= GrAphiMin, gphib= GrBphiMin), method="Newton") 
	if(echo)
		print(res$par)
	c(res$par, res$value )
}

respt <- t( sapply(1:NROW(initpt), NewLnsrch ) )

totalres <- cbind(1:NROW(initpt), initpt[, 1:2], NA, round(respt, 3))

#remove non convergent optimization sequences
finalres <- totalres[totalres[,9] == 0, -9]
colnames(finalres) <- c("num", "x1 init", "x2 init", "", "x1 final", "x2 final", "lam 1", "lam 2")

print(finalres)
print(trueGNE)

?GNE

#-------------------------------------------------------------------------------
#benchmark
z0 <- c(10, 10, 1, 1)

#min function
resMin <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMin), argjac=list(gphia= GrAphiMin, gphib= GrBphiMin), echo=FALSE)

resMin$compres

#FB function
resFB <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), echo=FALSE)

resFB$compres


#Mangasarian function
resMan <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMan, f=function(t) t^3), argjac=list(gphia= GrAphiMan, gphib= GrBphiMan, fprime=function(t) 3*t^2), echo=FALSE, control=list(maxit=200))

resMan$compres

#LT function
resLT <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiLT, q=4), argjac=list(gphia= GrAphiLT, gphib= GrBphiLT, q=4))

resLT$compres


#KK function
resKK <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiKK, lambda=3/2), argjac=list(gphia= GrAphiKK, gphib= GrBphiKK, lambda=3/2))

resKK$compres

	

}
\keyword{nonlinear}
\keyword{optimize}

