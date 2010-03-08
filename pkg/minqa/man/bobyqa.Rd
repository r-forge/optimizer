\name{bobyqa}
\alias{bobyqa}
\encoding{UTF-8}
\title{An R interface to the bobyqa implementation of Powell}
\description{
  The purpose of \code{bobyqa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
  Box constraints (bounds) on the parameters are permitted.
}
\usage{
bobyqa(par, fn, lower, upper, control = bobyqa.control(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates of the parameters
  of the objective function.  }
 \item{lower}{A numeric vector of lower bounds on the parameters. If the
 length is 1 the single lower bound is applied to all parameters. }
 \item{upper}{A numeric vector of upper bounds on the parameters.  If the
 length is 1 the single upper bound is applied to all parameters.}
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
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fval}{The value of the objective at the best set of parameters found.}
  \item{feval}{The number of function evaluations used.}
}
\references{
	M. J. D. Powell (2007)
	"Developments of NEWUOA for unconstrained minimization without
	derivatives",
	Cambridge University, Department of Applied Mathematics and
	Theoretical Physics,
	Numerical Analysis Group, Report NA2007/05,
        \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2007_05.pdf}.
	
	M. J. D. Powell (2009), "The BOBYQA algorithm for bound constrained
	optimization without derivatives", Report No. DAMTP 2009/NA06,
	Centre for Mathematical Sciences, University of Cambridge, UK.	
	
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
x1 <- bobyqa(fn=fr,par=c(1,2), lower=c(0,0), upper=c(4,4))
}
\keyword{nonlinear}
\keyword{optimize}

