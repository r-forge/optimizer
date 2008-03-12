\name{nlsolve}
\alias{nlsolve}
\title{Nonlinear System Solution by Minimization }
\description{Solve a nonlinear system by converting
  	  into a minimization problem and using optim's L-BFGS-B algorithm}
\usage{
     nlsolve(par, fn, lower=-Inf, upper=Inf, nstarts=10, ...)
 }
\arguments{
    \item{fn}{a function with a scalar real result (see details).}
    \item{par}{DOCUMENTa real scalar or vector argument to func, indicating the 
      point(s) at which the dfsaneient is to be calculated.}
    \item{lower}{DOCUMENT} 
    \item{upper}{DOCUMENT} 
    \item{nstarts}{DOCUMENT} 
    \item{...}{an additional arguments passed to \code{func}.}
}
\value{FOR EXAMPLEA real scalar or vector of the approximated nlsolveient(s).}

\details{ A solver for finding a zero of a system of non-linear equations
  Actually minimizes the squared-norm of the set of functions by calling optim()
  Uses the L-BFGS-B algorithm within optim()
  All the control parameters can be passed as in the call to optim()

  Author:  Ravi Varadhan, Center on Aging and Health, Johns Hopkins University, rvaradhan@jhmi.edu
  June 21, 2007

FOR EXAMPLE
   The function \code{nlsolve} calculates a numerical approximation of the 
   first derivative of \code{func} at the point \code{x}. Any additional 
   arguments in \dots are also passed to \code{func}, but the nlsolveient is not
   calculated with respect to these additional arguments.
   It is assumed \code{func}  is a scalar value function. If a vector \code{x} 
   produces a scalar
   result then \code{nlsolve} returns the numerical approximation of the nlsolveient
   at the point \code{x} (which has the same length as \code{x}).
   If a vector \code{x} produces a vector result then the result must have the
   same length as \code{x}, and it is assumed that this corresponds to applying
   the function to each of its arguments (for example, \code{sin(x)}). 
   In this case \code{nlsolve} returns the
   nlsolveient at each of the points in \code{x} (which also has the same length 
   as \code{x} - so be careful). An alternative for vector valued functions is
   provided by .
  
   If method is "simple", the calculation is done using a simple epsilon
   difference. For this case, only the element \code{eps} of \code{methods.args}
   is used.

   If method is "Richardson", the calculation
   is done by Richardson's extrapolation (see e.g. Linfield and Penny, 1989,
   or Fornberg and Sloan, 1994.)
   This method should be used if accuracy, as opposed to speed, is important.
   For this case, 
   \code{methods.args=list(eps=1e-4, d=0.01, r=4, show.details=FALSE)} 
   are used.
    \code{d} gives the fraction of \code{x} to use for the initial numerical 
      approximation. The default means the initial approximation uses
      \code{0.0001 * x}.
    \code{eps} is used instead of \code{d} for elements of \code{x} which are 
       zero.
    \code{r} gives the number of Richardson improvement iterations (repetitions
       with successly smaller \code{d}. The default \code{4} general provides 
       good results, but this can be increased to \code{6} for improved
       accuracy at the cost of more evaluations.
    \code{v} gives the reduction factor.
    \code{show.details} is a logical indicating if detailed calculations should 
        be shown.
   
  The general approach in the Richardson method is to iterate for \code{r} 
  iterations from initial 
  values for interval value \code{d},  using reduced factor \code{v}.
  The the first order approximation to the derivative with respect 
  to \eqn{x_{i}}{x_{i}} is

      \deqn{f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}{%
	    f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}
       
  This is repeated \code{r} times  with successively smaller \code{d}  and 
	  then Richardson extraplolation. is applied.

}
\references{ 
   FOR EXAMPLE
   Linfield, G. R. and Penny, J. E. T. (1989) \emph{Microcomputers in Numerical 
   Analysis}. New York: Halsted Press.

  Fornberg and Sloan (Acta Numerica, 1994, p. 203-267)
}
\seealso{
  \code{\link{spg}},
  \code{\link{dfsane}},
  \code{\link{sane}},
  \code{\link[stats]{numericDeriv}}
  \code{\link[stats]{optim}}
   }
\examples{
   expo1 <- function(x) {
    	n <- length(x)
    	f <- rep(NA, n)
    	f[1] <- exp(x[1] - 1) - 1
    	f[2:n] <- (2:n) * (exp(x[2:n] - 1) - x[2:n])
    	f
    	}

   nlsolve(par=runif(100), fn=expo1)
}
\keyword{multivariate}
