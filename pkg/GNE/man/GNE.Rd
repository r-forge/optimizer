\name{GNE}
\alias{GNE}
\title{GNE package}
\description{
		Generalized Nash Equilibrium computational methods.
	}
\usage{
GNE(approach = 
	c("non smooth", "fixed point", "minimization", "constrained equation"), 
	method = "default", xinit, control=list(), ...)

}
\arguments{
  \item{approach}{a character string for the approach: either \code{"non smooth"}, 
    \code{"fixed point"}, \code{"minimization"} or \code{"constrained equation"}.}
  \item{method}{a character string for the computation method: either \code{"default"} or the name
	of the method.}
  \item{xinit}{a numeric vector for the initial point.}
  \item{\dots}{further arguments to be passed to \code{GNE.nseq}, \code{GNE.fpeq} or \code{GNE.minpb}.}
  \item{control}{a list with control parameters.}
}	
\details{
Computing generalized Nash Equilibrium can be done in three different approaches.
\describe{
	\item{(i) extended KKT system}{It consists in solving the non smooth extended Karush-Kuhn-Tucker 
		(KKT) system \eqn{\Phi(z)=0}.}
	\item{(ii) fixed point approach}{It consists in solving equation \eqn{y(x)=x}.}
	\item{(iii) gap function minimization}{It consists in minimizing a gap function \eqn{min V(x)}.}
	\item{(iv) constrained equation}{It consists in solving \eqn{F(x)} such that \eqn{x}
	belongs to a specific set.}
}

The \code{GNE} function is a global function calling the appropriate function \code{\link{GNE.nseq}}, 
\code{\link{GNE.fpeq}}, \code{\link{GNE.ceq}} or \code{\link{GNE.minpb}}.
Benchmark functions comparing all methods for a given reformulation are 
available: see \code{\link{bench.GNE}}.


Additionnal utitilty functions are also available:
\code{\link{rejection}}, \code{\link{projector}}, \code{\link{stepfunc}},
 \code{\link{complementarity}} and \code{\link{funSSR}}.


TODO:
\describe{
\item{-}{write a pdf vignette.}
}

}
\value{
A list with components:
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
		\item{\code{fvec}}{a vector with function values.}
		\item{\code{approach}}{the name of the approach.}	
	}
}
\references{

 F. Facchinei, A. Fischer & V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.

 A. von Heusinger (2009),
 \emph{Numerical Methods for the Solution of the Generalized Nash Equilibrium Problem},
 Ph. D. Thesis.
 
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .


}
\seealso{

See \code{\link{GNE.fpeq}}, \code{\link{GNE.minpb}}, \code{\link{GNE.ceq}} 
and \code{\link{GNE.nseq}} for other approaches.

}
\author{
 Christophe Dutang
}
\keyword{math}


