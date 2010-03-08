\name{uobyqa}
\alias{uobyqa}
\encoding{UTF-8}
\title{An R interface to the uobyqa implementation of Powell}
\description{
  The purpose of \code{uobyqa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
 
}
\usage{
uobyqa(par, fn, control = uobyqa.control(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates. }
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{control}{
    An optional list of control settings.  See \code{\link{uobyqa.control}} for
    the names of the settable control values and their effect.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;  see \code{\link{uobyqa.control}} for
  details.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
}
\references{
     M. J. D. Powell, "The uobyqa software for unconstrained
     optimization 
     without derivatives",
     in \bold{Large-Scale Nonlinear Optimization},
     Series: Nonconvex Optimization and 
     Its Applications , Vol. 83, Di Pillo, Gianni; Roma, Massimo (Eds.)
     2006, New York: Springer US.

     M. J. D. Powell,
     "Developments of uobyqa for minimization without derivatives",
     IMA Journal of  Numerical Analysis, 2008; 28: 649-664.

     
     Description was taken from comments in the Fortran code of
     M. J. D. Powell on which \pkg{minqa} is based.
}
\seealso{\code{\link{optim}}, \code{\link{nlminb}}}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
x3 <- uobyqa(fn=fr,par=c(1,2))
}
\keyword{nonlinear}
\keyword{optimize}

