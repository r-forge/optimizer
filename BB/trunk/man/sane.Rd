\name{sane}
\alias{sane}
\title{Spectral-Gradient Solution of Nonlinear System}
\description{Spectral-gradient based solution of nonlinear systems}
\usage{
     sane(p, func, M=5, maxiter=1000, lower=-Inf, upper=Inf, tol=1.e-07, 
       trace=FALSE, ...)
 }
\arguments{
    \item{func}{a function with a scalar real result (see details).}
    \item{p}{DOCUMENTa real scalar or vector argument to func, indicating the 
      point(s) at which the dfsaneient is to be calculated.}
    \item{M}{DOCUMENT} 
    \item{maxiter}{DOCUMENT} 
    \item{lower}{DOCUMENT} 
    \item{upper}{DOCUMENT} 
    \item{tol}{DOCUMENT} 
    \item{trace}{DOCUMENT} 
    \item{...}{an additional arguments passed to \code{func}.}
}
\value{FOR EXAMPLEA real scalar or vector of the approximated saneient(s).}

\details{Non-monotone spectral method for finding a root of nonlinear systems

R translation (with minor modifications):  Ravi Varadhan, Johns Hopkins University, February 23, 2008.

FOR EXAMPLE
   The function \code{sane} calculates a numerical approximation of the 
   first derivative of \code{func} at the point \code{x}. Any additional 
   arguments in \dots are also passed to \code{func}, but the saneient is not
   calculated with respect to these additional arguments.
   It is assumed \code{func}  is a scalar value function. If a vector \code{x} 
   produces a scalar
   result then \code{sane} returns the numerical approximation of the saneient
   at the point \code{x} (which has the same length as \code{x}).
   If a vector \code{x} produces a vector result then the result must have the
   same length as \code{x}, and it is assumed that this corresponds to applying
   the function to each of its arguments (for example, \code{sin(x)}). 
   In this case \code{sane} returns the
   saneient at each of the points in \code{x} (which also has the same length 
   as \code{x} - so be careful). An alternative for vector valued functions is
   provided by.
  
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
  La Cruz and Raydan, Optim Methods and Software 2003, 18 (583-599)
}
\seealso{
  \code{\link{spg}},
  \code{\link{dfsane}},
  \code{\link{nlsolve}},
  \code{\link[stats]{numericDeriv}}
   }
\examples{
   expo1 <- function(x) {
    	n <- length(x)
    	f <- rep(NA, n)
    	f[1] <- exp(x[1] - 1) - 1
    	f[2:n] <- (2:n) * (exp(x[2:n] - 1) - x[2:n])
    	f
    	}

   sane(p=runif(100), func=expo1)
}
\keyword{multivariate}
