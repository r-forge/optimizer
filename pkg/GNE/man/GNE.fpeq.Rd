\name{GNE.fpeq}
\alias{GNE.fpeq}
\title{Fixed point equation reformulation of the GNE problem.}
\description{
		Fixed point equation reformulation via the NI function of the GNE problem.
	}
\usage{
fpeq(xinit,	yfunc, stepfunc, Vfunc, method=c("pure","UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"), 
	control=list(), ...)

}
\arguments{
  \item{xinit}{initial point.}
  \item{yfunc}{the y function.}
  \item{stepfunc}{the step function, only needed when \code{method="UR"}.}
  \item{Vfunc}{the value function used to compute the merit function.}
  \item{method}{either \code{"pure"}, \code{"UR"}, \code{"vH"}, \code{"RRE"}, \code{"MPE"}, 
	\code{"SqRRE"} or \code{"SqMPE"} method, see details.}
  \item{control}{a list with control parameters.}
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions.}
  
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
		\item{\code{value}}{The value of the merit function.}
		\item{\code{counts}}{A two-element integer vector giving the number of 
			calls to \code{yfunc} and \code{Vfunc}/\code{stepfunc} respectively.}					
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
 
 S. Uryasev & R.Y. Rubinstein (1994), 
 \emph{On relaxation algorithms in computation of noncooperative equilibria}, 
 IEEE Transactions on Automatic Control.
 
 R. Varadhan & C. Roland (2004),
 \emph{Squared Extrapolation Methods (SQUAREM): A New Class of Simple and Efficient Numerical 
 Schemes for Accelerating the Convergence of the EM Algorithm},
 Johns Hopkins University, Dept. of Biostatistics Working Papers.
 
}

\seealso{

See \code{\link{GNE.nseq}} and \code{\link{GNE.min}} 
for other approaches.


}
\author{
 Christophe Dutang
}
\keyword{nonlinear}
\keyword{optimize}

