\name{uobyqa}
\alias{uobyqa}
\encoding{UTF-8}
\title{An R interface to the uobyqa implementation of Powell}
\description{
  The purpose of \code{uobyqa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
  There can be some freedom in the interpolation conditions, which is
  taken up by minimizing the Frobenius norm of the change to the second
  derivative of the quadratic model, beginning with a zero matrix.
  (Description taken from comments in Fortran code of M J D Powell.)
}
\usage{
uobyqa(par, fn, control = uobyqa.control(), check.hessian = TRUE, \dots)
}
\arguments{
 \item{par}{A list or numeric vector of starting estimates. If
   \code{par} is a list, then each element must be of length 1. }
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{control}{
    An optional list of control settings.  See \code{\link{uobyqa.control}} for
    the names of the settable control values and their effect.
 }
 \item{check.hessian}{Set TRUE to check if Hessian is positive definite. Sundar:
   Is this an addition you made? It does not appear in the original Fortran.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;  see \code{\link{uobyqa.control}} for
    details.
    Successful completion.
    \code{uobyqa} for R is still a work in progress, so users should watch
    the console output.
    Unsuccessful completion.
    Unsuccessful termination of \code{nls.lm} can be due to improper
    input parameters, arithmetic interrupts, an excessive number of
    function evaluations, or an excessive number of iterations. 
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
}
\references{
     M.J.D. Powell, "The uobyqa software for unconstrained optimization without derivatives",
     in \bold{Large-Scale Nonlinear Optimization}, Series: Nonconvex Optimization and 
     Its Applications , Vol. 83, Di Pillo, Gianni; Roma, Massimo (Eds.)
     2006, New York: Springer US

     M.J.D. Powell, "Developments of uobyqa for minimization without derivatives"
     IMA J Numer Anal.2008; 28: 649-664 
}
\seealso{\code{\link{optim}}}
\examples{
# There are no examples yet.
}
\keyword{nonlinear}
\keyword{optimize}

