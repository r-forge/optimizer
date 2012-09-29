\name{hessian}
\alias{hessian}
\alias{hessian.default}

\title{Calculate Hessian Matrix}
\description{Calculate a numerical approximation to the Hessian matrix of a 
function at a parameter value.}
\usage{
    hessian(func, x, method="Richardson", method.args=list(), ...)

    \method{hessian}{default}(func, x, method="Richardson",
        method.args=list(), ...)
}
\arguments{
    \item{func}{a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{x}{the parameter vector first argument to func.}
    \item{method}{one of \code{"Richardson"} or \code{"complex"} indicating 
    the method to use for the approximation.}
    \item{method.args}{arguments passed to method.  See \code{\link{grad}}. 
       (Arguments not specified remain with their default values.)}
    \item{...}{an additional arguments passed to \code{func}.
          WARNING: None of these should have names matching other arguments of this function.}
}
\value{An n by n matrix of the Hessian of the function calculated at the 
       point \code{x}.}
\details{
  The function \code{hessian} calculates an numerical approximation to 
  the n x n second derivative of a scalar real valued function with n-vector
  argument. 
  
  The argument \code{method} can be \code{"Richardson"} or \code{"complex"}.
  Method \code{"simple"} is not supported. Methods \code{"complex"} is
  described more in \code{\link{grad}}.
  For method "complex" \code{methods.args=list(eps=.Machine$double.eps)} is the
  default and only \code{eps} is used by this method.

  Methods \code{"Richardson"} is uses \code{\link{genD}} and extracts the 
  second derivative. For this method 
   \code{methods.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), 
   r=4, v=2, show.details=FALSE)} is set as the default.
}
\seealso{
    \code{\link{jacobian}},
    \code{\link{grad}},
    \code{\link{genD}}
}
\keyword{multivariate}
