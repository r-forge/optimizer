\name{GNE.min}
\alias{GNE.min}
\title{Non smooth equation reformulation of the GNE problem.}
\description{
		Non smooth equation reformulation via the extended KKT system of the GNE problem.
	}
\usage{
GNE.min(xinit, gap, gradgap, arggap=list(), arggrad=list(), method, control=list(), ...)
}
\arguments{
  \item{xinit}{initial point.}
  \item{gap}{the gap/value function.}
  \item{gradgap}{the gradient of the gap/value function.}
  \item{arggap}{further arguments to be passed to \code{phi}.}
  \item{arggrad}{further arguments to be passed to \code{jacphi}.}
  \item{method}{either \code{"BB"} or \code{BFGS}.}
  \item{control}{a list with control parameters.}
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions \code{phi} and \code{jacphi}.}
  
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
		\item{\code{value}}{The value of the merit function.}
		\item{\code{outer.counts}}{A two-element integer vector giving the number of 
			calls to \code{Gap} and \code{gradGap} respectively.}				
		\item{\code{outer.iter}}{The outer iteration number.}
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
		\item{\code{inner.iter.fn}, \code{inner.iter.gr}}{The iteration number
			for the function and the gradient when computing the gap function or
			its gradient (if appropriate).}	
		\item{\code{inner.counts.fn}, \code{inner.counts.gr}}{A two-element integer 
			vector giving the number of calls to the function and the gradient 
			when computing the gap function or its gradient (if appropriate).}			
		\item{\code{message}}{a string describing the termination code}	
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

See \code{\link{GNE.fpeq}} and \code{\link{GNE.nseq}} 
for other approaches.

}
\author{
 Christophe Dutang
}
\keyword{nonlinear}
\keyword{optimize}

