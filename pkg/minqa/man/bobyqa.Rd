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
bobyqa(par, fn, lower = -Inf, upper = Inf, control = list(), \dots)
}
\arguments{
  \item{par}{A numeric vector of starting estimates of the parameters
    of the objective function.}
  \item{fn}{A function that returns the value of the objective at the
    supplied set of parameters \code{par} using auxiliary data in \dots.
    The first argument of \code{fn} must be \code{par}. }
  \item{lower}{A numeric vector of lower bounds on the parameters. If the
    length is 1 the single lower bound is applied to all parameters.}
  \item{upper}{A numeric vector of upper bounds on the parameters.  If the
    length is 1 the single upper bound is applied to all parameters.}
  \item{control}{
    An optional list of control settings.  See the details section for
    the names of the settable control values and their effect.}
  \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  The function \code{fn} must return a scalar numeric value.

  The \code{control} argument is a list.  Possible named values in the
  list and their defaults are:
  \describe{
    \item{npt}{
      The number of points used to approximate the objective function
      via a quadratic approximation. The value of npt must be in the
      interval \eqn{[n+2,(n+1)(n+2)/2]} where \eqn{n} is the number of
      parameters in \code{par}. Choices that exceed \eqn{2*n+1} are not
      recommended.  If not defined, it will be set to \eqn{\min(n * 2,
	n+2)}{min(n * 2, n+2)}.
    }
    \item{rhobeg}{
      \code{rhobeg} and \code{rhoend} must be set to the initial and final
      values of a trust region radius, so both must be positive with
      \code{0 < rhoend < rhobeg}. Typically \code{rhobeg} should be about
      one tenth of the greatest expected change to a variable.  If the
      user does not provide a value, this will be set to the maximum box
      constraint interval (for finite limits), i.e., \code{max(xu-xl)} or
      to \code{max(par) / 2)} if the bounds are infinite (i.e., no box
      constraints are specified).
    }
    \item{rhoend}{
      The smallest value of the trust region radius that is allowed. If
      not defined, then 1e-6 times the value set for \code{rhobeg} will be
      used.
    }
    \item{iprint}{
      The value of \code{iprint} should be set to \code{0, 1, 2 or 3},
      which controls the amount of printing. Specifically, there is no
      output if \code{iprint=0} and there is output only at the return if
      \code{iprint=1}. Otherwise, each new value of \code{rho} is printed,
      with the best vector of variables so far and the corresponding value
      of the objective function. Further, each new value of the objective
      function with its variables are output if \code{iprint=3}.  Default
      value is \code{0}.
    }
    \item{maxfun}{
      The maximum allowed number of function evaluations. If this is
      exceeded, the method will terminate.
    } 
  }
    
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
  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf}.
  
  Description was taken from comments in the Fortran code of
  M. J. D. Powell on which \pkg{minqa} is based.
}
\seealso{\code{\link{optim}}, \code{\link{nlminb}}}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
}
(x1 <- bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4)))
## => optimum at c(1, 1) with fval = 0
}
\keyword{nonlinear}
\keyword{optimize}

