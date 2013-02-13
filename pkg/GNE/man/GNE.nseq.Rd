\name{GNE.nseq}
\alias{GNE.nseq}
\alias{bench.GNE.nseq}
\alias{bench.GNE.nseq.LM}
\title{Non smooth equation reformulation of the GNE problem.}
\description{
		Non smooth equation reformulation via the extended KKT system of the GNE problem.
	}
\usage{
GNE.nseq(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="default", control=list(), silent=TRUE, ...)
	
#bench.GNE.nseq(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), echo=FALSE, ...)
#bench.GNE.nseq.LM(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), echo=FALSE, control=list())
}
\arguments{
  \item{init}{Initial values for the parameters to be optimized over: \eqn{z=(x, lambda, mu)}.}
  \item{dimx}{a vector of dimension for \eqn{x}.}
  \item{dimlam}{a vector of dimension for \eqn{lambda}.}
  \item{grobj}{gradient of the objective function (to be minimized), see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}
  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}	

  \item{constr}{constraint function (\eqn{g^i(x)<=0}), see details.}
  \item{argconstr}{a list of additional arguments of the constraint function.}
  \item{grconstr}{gradient of the constraint function, see details.}
  \item{arggrconstr}{a list of additional arguments of the constraint gradient.}
  \item{heconstr}{Hessian of the constraint function, see details.}
  \item{argheconstr}{a list of additional arguments of the constraint Hessian.}

  \item{compl}{the complementarity function with (at least) two arguments: \code{compl(a,b)}.}
  \item{argcompl}{list of possible additional arguments for \code{compl}.}
  \item{gcompla}{derivative of the complementarity function w.r.t. the first argument.}
  \item{gcomplb}{derivative of the complementarity function w.r.t. the second argument.}   
  
  \item{dimmu}{a vector of dimension for \eqn{mu}.}
  \item{joint}{joint function (\eqn{h(x)<=0}), see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{grjoint}{gradient of the joint function, see details.}
  \item{arggrjoint}{a list of additional arguments of the joint gradient.}
  \item{hejoint}{Hessian of the joint function, see details.}
  \item{arghejoint}{a list of additional arguments of the joint Hessian.}

  \item{method}{a character string specifying the method \code{"Newton"},
	\code{"Broyden"}, \code{"Levenberg-Marquardt"} or \code{"default"} 
	which is \code{"Newton"}. }
  \item{control}{a list with control parameters.}
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions \code{phi} and \code{jacphi}.}
  \item{silent}{a logical to get some traces. Default to \code{FALSE}.}	
}
\details{

Functions in argument must respect the following template:
\itemize{
  \item{\code{constr} must have arguments the current iterate \code{z}, the player number \code{i} 
  and optionnally additional arguments given in a list.}
  \item{\code{grobj}, \code{grconstr} must have arguments the current iterate \code{z}, the player number \code{i}, 
  the derivative index \code{j} and optionnally additional arguments given in a list.}
  \item{\code{heobj}, \code{heconstr} must have arguments the current iterate \code{z}, the player number \code{i}, 
  the derivative indexes \code{j}, \code{k} and optionnally additional arguments given in a list.}
  \item{\code{compl}, \code{gcompla}, \code{gcomplb} must have two arguments \code{a}, \code{b} and optionnally additional arguments given in a list.}
  \item{\code{joint} must have arguments the current iterate \code{z}
  and optionnally additional arguments given in a list.}
  \item{\code{grjoint} must have arguments the current iterate \code{z},
  the derivative index \code{j} and optionnally additional arguments given in a list.}
  \item{\code{hejoint} must have arguments the current iterate \code{z}, 
  the derivative indexes \code{j}, \code{k} and optionnally additional arguments given in a list.}
}

\code{GNE.nseq} solves the GNE problem via a non smooth reformulation of the KKT system. 
\code{bench.GNE.nseq} carries out a benchmark of the computation methods (Newton and Broyden 
direction with all possible global schemes) for a given initial point.
\code{bench.GNE.nseq.LM} carries out a benchmark of the Levenberg-Marquardt computation method.

This approach consists in solving the extended Karush-Kuhn-Tucker 
(KKT) system denoted by \eqn{\Phi(z)=0}{Phi(z)=0}, where eqn{z} is formed by the players strategy 
\eqn{x} and the Lagrange multiplier \eqn{\lambda}{lambda}.
The root problem \eqn{\Phi(z)=0}{Phi(z)=0} is solved by an iterative scheme \eqn{z_{n+1} = z_n + d_n},
where the direction \eqn{d_n} is computed in three different ways. Let \eqn{J(x)=Jac\Phi(x)}{J(x)=Jac Phi(x)}.
\describe{
	\item{(a) Newton:}{The direction solves the system \eqn{J(z_n) d = - \Phi(z_n)
		}{J(z_n) d = - Phi(z_n)}, 
		generally called the Newton equation.}
	\item{(b) Broyden:}{It is a quasi-Newton method aiming to solve an approximate version
		of the Newton equation \eqn{d = -\Phi(z_n) W_n}{d = -Phi(z_n) W_n} where \eqn{W_n} is computed
		by an iterative scheme. In the current implementation, \eqn{W_n} is updated
		by the Broyden method. }
	\item{(c) Levenberg-Marquardt:}{The direction solves the system 
		\deqn{ \left[ J(z_n)^T J(z_n) + \lambda_n^\delta I \right] d = - J(z_n)^T\Phi(x_n)
			}{ [J(z_n)^T J(z_n) + lambda_n^delta I] d = - J(z_n)^T Phi(z_n),}
		where \eqn{I} denotes the identity matrix, \eqn{\delta}{delta} is a parameter in [1,2] 
		and \eqn{\lambda_n = ||\Phi(z_n)|| }{lambda_n = ||Phi(z_n)||} if \code{LM.param="merit"},
		\eqn{||J(z_n)^T \Phi(z_n)|| }{||J(z_n)^T Phi(z_n)||} if \code{LM.param="jacmerit"},
		the minimum of both preceding quantities if \code{LM.param="min"}, or an adatpive
		parameter according to Fan(2003) if \code{LM.param="adaptive"}.
		}
}

In addition to the computation method, a globalization scheme can be choosed using the \code{global}
argument, via the \code{\dots} argument. Available schemes are 
\describe{
	\item{(1) Line search:}{ if \code{global} is set to \code{"qline"} or \code{"gline"}, a line search
		is used with the merit function being half of the L2 norm of \eqn{Phi}, respectively with a
		quadratic or a geometric implementation.}
	\item{(2) Trust region:}{ if \code{global} is set to \code{"dbldog"} or \code{"pwldog"}, a trust
		region is used respectively with a double dogleg or a Powell (simple) dogleg implementation.
		This global scheme is not available for the Levenberg-Marquardt direction.}
	\item{(3) None:}{ if  \code{global} is set to \code{"none"}, no globalization is done. }
}	
The default value of \code{global} is \code{"gline"}. Note that in the special case of 
the Levenberg-Marquardt direction with adaptive parameter, the global scheme must be \code{"none"}.

In the GNEP context, details on the methods can be found in Facchinei, Fischer & Piccialli (2009), \code{"Newton"}
corresponds to method 1 and \code{"Levenberg-Marquardt"} to method 3. In a general nonlinear
equation framework, see Dennis & Moree (1977), Dennis & Schnabel (1996) or Nocedal & Wright (2006), 


The implementation relies heavily on the 
\code{\link{nleqslv}} function of the package of the same name. So full details on the control parameters are
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
	\item{\code{delta}}{A numeric \code{delta} in [1, 2], default to 2, for
	 the Levenberg-Marquardt method only.
	}
	\item{\code{LM.param}}{A character string, default to \code{"merit"}, for
	 the Levenberg-Marquardt method only.
	}
	\item{\code{maxit}}{The maximum number of major iterations. The default value is 150 if a 
	global strategy has been specified.}
	\item{\code{trace}}{Non-negative integer. A value of 1 will give a detailed report of the
	progress of the iteration, default 0.}
}

\code{\dots} are further arguments to be passed to the optimization routine, 
that is \code{global}, \code{xscalm}, \code{silent}. See above for the globalization scheme. 
The \code{xscalm} is a scaling parameter to used, either \code{"fixed"} (default)
or \code{"auto"}, for which scaling factors are calculated from the euclidean norms of the 
columns of the jacobian matrix. See \code{\link{nleqslv}} for details.
The \code{silent} argument is a logical to report or not the optimization process, default
to \code{FALSE}.
 

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
		\item{\code{message}}{a string describing the termination code.}	
		\item{\code{fvec}}{a vector with function values.}
	}
\code{bench.GNE.nseq} returns a list with components:
	\describe{
		\item{\code{compres}}{a data.frame summarizing the different computations.}
		\item{\code{reslist}}{a list with the different results from \code{GNE.nseq}.}
	} 	

}
\references{

 J.E. Dennis and J.J. Moree (1977), 
 \emph{Quasi-Newton methods, Motivation and Theory},
 SIAM review.

 J.E. Dennis and R.B. Schnabel (1996), 
 \emph{Numerical methods for unconstrained optimization and nonlinear equations},
 SIAM.
 
 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
 
 J.-Y. Fan (2003),
 \emph{A modified Levenberg-Marquardt algorithm for singular 
 system of nonlinear equations},
 Journal of Computational Mathematics.
 
 B. Hasselman (2011), 
 \emph{nleqslv: Solve systems of non linear equations},
 R package.

 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem 
 using Nikaido-Isoda-type functions},
 Comput Optim Appl .

 J. Nocedal and S.J. Wright (2006), 
 \emph{Numerical Optimization}, 
 Springer Science+Business Media

}
\seealso{

See \code{\link{GNE.fpeq}} and \code{\link{GNE.min}} 
for other approaches; \code{\link{funSSR}} and
\code{\link{jacSSR}} for template functions of \eqn{\Phi} and \eqn{Jac\Phi} and
\code{\link{complementarity}} for complementarity functions.

See also \code{\link{nleqslv}} for some optimization details.
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



#true value is (3/4, 1/4, 1/2, 1/2)

z0 <- rep(0, sum(dimx)+sum(dimlam))

funSSR(z0, dimx, dimlam, grobj=grobj, constr=g, grconstr=grg, compl=phiFB, echo=FALSE)

	
jacSSR(z0, dimx, dimlam, heobj=heobj, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, NULL, heobj=heobj, NULL, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, NULL, heobj=heobj, NULL, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1))



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

z0 <- rep(0, sum(dimx)+sum(dimlam))

funSSR(z0, dimx, dimlam, grobj=grobj, myarg, constr=g, grconstr=grg, compl=phiFB, echo=FALSE)

jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1))


	

}
\keyword{nonlinear}
\keyword{optimize}

