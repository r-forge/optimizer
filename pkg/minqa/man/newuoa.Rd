\name{newuoa}
\alias{newuoa}
\encoding{UTF-8}
\title{An R interface to the NEWUOA implementation of Powell}
\description{
  The purpose of \code{newuoa} is to minimize a function of many variables 
  by a trust region method that forms quadratic models by interpolation.
 
}
\usage{
newuoa(par, fn, control = list(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates.  }
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{control}{
    An optional list of control settings.  See the details section for
    the names of the settable control values and their effect.
}
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;

  Powell's Fortran code has been slightly modified (thanks to Doug Bates for 
  help on this) to avoid use of PRINT statements. Output is now via calls to
  C routines set up to work with the routines BOBYQA, NEWUOA and UOBYQA.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fval}{The value of the objective at the best set of parameters found.}
  \item{feval}{Number of function evaluations to determine the optimum}
}
\references{
  M. J. D. Powell, "The NEWUOA software for unconstrained optimization
  without derivatives",
  in \bold{Large-Scale Nonlinear Optimization}, Series: Nonconvex
  Optimization and Its Applications , Vol. 83, Di Pillo, Gianni;
  Roma, Massimo (Eds.) 2006, New York: Springer US.

  M. J. D. Powell,
  "Developments of NEWUOA for minimization without derivatives"
  IMA Journal of Numerical Analysis, 2008; 28: 649-664. 
  
  M. J. D. Powell (2007)
  "Developments of NEWUOA for unconstrained minimization without
  derivatives"
  Cambridge University, Department of Applied Mathematics and
  Theoretical Physics,
  Numerical Analysis Group, Report NA2007/05,
  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2007_05.pdf}.

  Description was taken from comments in the Fortran code of
  M. J. D. Powell on which \pkg{minqa} is based.
}
\seealso{\code{\link{optim}}, \code{\link{nlminb}}}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
}
(x2 <- newuoa(c(1, 2), fr))
## => optimum at c(1, 1) with fval = 0
}
\keyword{nonlinear}
\keyword{optimize}

