\name{GNE.ceq}
\alias{GNE.ceq}
\title{Constrained equation reformulation of the GNE problem.}
\description{
		Constrained equation reformulation via the extended KKT system of the GNE problem.
	}
\usage{
GNE.ceq(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="IP", control=list(), silent=TRUE, ...)

}
\arguments{
  \item{init}{Initial values for the parameters to be optimized over: \eqn{z=(x, lambda, mu)}.}
  \item{dimx}{a vector of dimension for \eqn{x}.}
  \item{dimlam}{a vector of dimension for \eqn{lambda}.}
  
  \item{grobj}{gradient of the objective function, see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}
  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}	
  
  \item{constr}{constraint function, see details.}
  \item{argconstr}{a list of additional arguments of the constraint function.}
  \item{grconstr}{gradient of the constraint function, see details.}
  \item{arggrconstr}{a list of additional arguments of the constraint gradient.}
  \item{heconstr}{Hessian of the constraint function, see details.}
  \item{argheconstr}{a list of additional arguments of the constraint Hessian.}

  \item{dimmu}{a vector of dimension for \eqn{mu}.}
  \item{joint}{joint function, see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{grjoint}{gradient of the joint function, see details.}
  \item{arggrjoint}{a list of additional arguments of the joint gradient.}
  \item{hejoint}{Hessian of the joint function, see details.}
  \item{arghejoint}{a list of additional arguments of the joint Hessian.}

  \item{method}{a character string specifying the method \code{"Newton"},
	\code{"Broyden"}, \code{"Levenberg-Marquardt"}. }
  \item{control}{a list with control parameters.}
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions \code{H} and \code{jacH}.}
  \item{silent}{a logical to get some traces. Default to \code{FALSE}.}	
}
\details{
\code{GNE.ceq} solves the GNE problem via a constrained equation reformulation of the KKT system. 

This approach consists in solving the extended Karush-Kuhn-Tucker 
(KKT) system denoted by \eqn{H(z)=0}, where eqn{z} is formed by the players strategy 
\eqn{x}, the Lagrange multiplier \eqn{\lambda}{lambda} and the slate variable \eqn{w}.
The root problem \eqn{H(z)=0} is solved by an iterative scheme \eqn{z_{n+1} = z_n + d_n},
where the direction \eqn{d_n} is computed in two different ways. Let \eqn{J(x)=Jac H(x)}.
\describe{
	\item{(a) modified Newton:}{The direction solves the system \eqn{H(z_n) + J(x_n) d = sigma_n a^T H(z_n) / ||a||_2^2 a
		}.}
	\item{(b) ?:}{}
}

In addition to the computation method, a globalization scheme can be choosed using the \code{global}
argument, via the \code{\dots} argument. Available schemes are 
\describe{
	\item{(1) Line search:}{ if \code{global} is set to \code{"qline"} or \code{"gline"}, a line search
		is used with the merit function being half of the L2 norm of \eqn{Phi}, respectively with a
		quadratic or a geometric implementation.}
	\item{(2) None:}{ if  \code{global} is set to \code{"none"}, no globalization is done. }
}	
The default value of \code{global} is \code{"gline"}. 

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
	\item{\code{btol}}{The backtracking tolerance.
	The	default value is \eqn{10^{-2}}{1e-2}.
	}	
	\item{\code{maxit}}{The maximum number of major iterations. The default value is 100 if a 
	global strategy has been specified.}
	\item{\code{trace}}{Non-negative integer. A value of 1 will give a detailed report of the
	progress of the iteration, default 0.}
	\item{\code{sigma}, \code{delta}, \code{zeta}}{Parameters initialized to \code{1/2}, 
	\code{1}, \code{length(init)/2}.}
	\item{\code{forcingpar}}{Forcing parameter set to 0.1.}
}

}
\value{
\code{GNE.ceq} returns a list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{value}}{The value of the merit function.}
		\item{\code{counts}}{A two-element integer vector giving the number of calls to 
			\code{H} and \code{jacH} respectively.}			
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
		\item{\code{message}}{a string describing the termination code.}	
		\item{\code{fvec}}{a vector with function values.}
	}

}
\references{

 J.E. Dennis and J.J. Moree (1977), 
 \emph{Quasi-Newton methods, Motivation and Theory},
 SIAM review.

 A. Dreves, F. Facchinei, C. Kanzow and S. Sagratella (2011), 
 \emph{On the solutions of the KKT conditions of generalized Nash equilibrium problems}, 
 SIAM Journal on Optimization 21(3), 1082–1108.


}
\seealso{

See \code{\link{GNE.fpeq}} and \code{\link{GNE.min}} 
for other approaches; \code{\link{funCER}} and
\code{\link{jacCER}} for template functions of \eqn{H} and \eqn{Jac H}.

}
\author{
 Christophe Dutang
}
\examples{


#-------------------------------------------------------------------------------
# (1) Example 5 of von Facchinei et al. (2007)
#-------------------------------------------------------------------------------

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j)
{
	if(i == 1)
		res <- c(2*(x[1]-1), 0)
	if(i == 2)
		res <- c(0, 2*(x[2]-1/2))
	res[j]	
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k)
	2 * (i == j && j == k)

dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	sum(x[1:2]) - 1
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	1
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0


x0 <- rep(0, sum(dimx))
z0 <- c(x0, 2, 2, max(10, 5-g(x0, 1) ), max(10, 5-g(x0, 2) ) )

#true value is (3/4, 1/4, 1/2, 1/2)
GNE.ceq(z0, dimx, dimlam, grobj=grobj, heobj=heobj, 
	constr=g, grconstr=grg, heconstr=heg, method="IP", 
	control=list(trace=0, maxit=10))




#-------------------------------------------------------------------------------
# (2) Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------


#constants
myarg <- list(d= 20, lambda= 4, rho= 1)

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	res <- -arg$rho * x[i]
	if(i == j)
	res <- res + arg$d - arg$lambda - arg$rho*(x[1]+x[2])
	-res
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
	arg$rho * (i == j) + arg$rho * (j == k)	


dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	-x[i]
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	-1*(i == j)
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0

#true value is (16/3, 16/3, 0, 0) 

x0 <- rep(0, sum(dimx))
z0 <- c(x0, 2, 2, max(10, 5-g(x0, 1) ), max(10, 5-g(x0, 2) ) )


GNE.ceq(z0, dimx, dimlam, grobj=grobj, heobj=heobj, arggrobj=myarg, 
	argheobj=myarg, constr=g, grconstr=grg, heconstr=heg,
	method="IP", control=list(trace=0, maxit=10))


	

}
\keyword{nonlinear}
\keyword{optimize}

