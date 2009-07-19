\name{bobyqa}
\alias{bobyqa}
\encoding{UTF-8}
\title{An R interface to the bobyqa implementation of Powell}
\description{
  The purpose of \code{bobyqa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
  Box constraints (bounds) on the parameters are permitted.
  There can be some freedom in the interpolation conditions, which is
  taken up by minimizing the Frobenius norm of the change to the second
  derivative of the quadratic model, beginning with a zero matrix.
  (Description taken from comments in Fortran code of M J D Powell.)
}
\usage{
bobyqa(par, xl, xu, fn, control = bobyqa.control(), \dots)
}
\arguments{
 \item{par}{A list or numeric vector of starting estimates of the parameters
  of the objective function. If \code{par} is a list, then each element 
  must be of length 1. }
 \item{xl}{A list or numeric vector of lower bounds on the parameters.}
 \item{xu}{A list or numeric vector of upper bounds on the parameters.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{control}{
    An optional list of control settings.  See \code{\link{bobyqa.control}} for
    the names of the settable control values and their effect.
}
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;  see \code{\link{bobyqa.control}} for
    details.
    Successful completion.
    \code{bobyqa} for R is still a work in progress, so users should watch
    the console output.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
}
\references{
-------
}
\seealso{\code{\link{optim}}}
\examples{
# There are no examples yet.
}
\keyword{nonlinear}
\keyword{optimize}

