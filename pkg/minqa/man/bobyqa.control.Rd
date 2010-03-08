\name{bobyqa.control}
\Rdversion{1.1}
\alias{bobyqa.control}
\title{
Controls for an R interface to the bobyqa implementation of Powell
}
\description{
  The purpose of \code{bobyqa.control} is to set control parameters for
  the \code{bobyqa} Fortran code of M. J. D. Powell. 
}
\usage{
bobyqa.control(npt = NA, rhobeg = NA, rhoend = NA, iprint = 0,
maxfun = 10000, obstop=TRUE, force.start=FALSE)
}

\arguments{
  \item{npt}{
   The number of points used to approximate the objective function via
   a quadratic approximation. The value of npt must be in the interval 
   \eqn{[n+2,(n+1)(n+2)/2]} where 
 \eqn{n}
 is the number of parameters in \code{par}. Choices that exceed
 \eqn{2*n+1} are 
 not recommended.
 If not defined (i.e., the default \code{NA} is used), it will be set to
 \eqn{\min(n * 2, n+2)}{min(n * 2, n+2)}. }
\item{rhobeg}{
   \code{rhobeg} and \code{rhoend}
    must be set to the initial and final values of a trust
    region radius, so both must be positive with
    \code{0 < rhoend < rhobeg}. Typically
    \code{rhobeg} should be about one tenth of the greatest expected
    change to a variable.
     If the user does not provide a value, this will be set
  to the maximum box constraint interval (for finite limits), i.e.,
        \code{max(xu-xl)}
   or to 
        \code{max(par) / 2)}
   if the bounds are infinite (i.e., no box constraints are specified).
}
  \item{rhoend}{
    The smallest value of the
    trust region radius that is allowed. If not defined, then 
   1e-6 times the value set for \code{rhobeg} will be used.
}
  \item{iprint}{
 The value of \code{iprint} should be set to \code{0, 1, 2 or 3},
  which controls the
  amount of printing. Specifically, there is no output if \code{iprint=0} and
       there is output only at the return if \code{iprint=1}. Otherwise, each
       new
       value of \code{rho}
       is printed, with the best vector of variables so far and
       the corresponding value of the objective function. Further, each new
       value of the objective function with its variables are output if
       \code{iprint=3}.
       Default value is \code{0}.
}
  \item{maxfun}{
  The maximum allowed number of function evaluations. If this is exceeded, the
  method will terminate.
}
 
  \item{obstop}{
  Set \code{TRUE} if out of bounds parameters force program to stop. Otherwise
	parameters are reset to nearest bound.
}
  \item{force.start}{
    Powell tries to move parameters at least \code{rhobeg} from the nearest
    bound. This is the default action (\code{force.start=FALSE)}. When
    \code{TRUE}, use user supplied start regardless.
      }
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;  see \code{\link{bobyqa.control}} for
    details.
}

\keyword{optimize}

