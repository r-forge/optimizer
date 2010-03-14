\name{newuoa.control}
\Rdversion{1.1}
\alias{newuoa.control}
\title{
Controls for an R interface to the NEWUOA implementation of Powell
}
\description{
  The purpose of \code{newuoa.control} is to set control parameters for
  the NEWUOA Fortran code of M. J. D. Powell. 
}
\usage{
newuoa.control(npt = NA, rhobeg = NA, rhoend = NA, iprint = 0,
maxfun = 10000)
}
\arguments{
  \item{npt}{
   The number of points used to approximate the objective function via
   a quadratic approximation. The value of npt must be in the interval 
   \eqn{[n+2,(n+1)(n+2)/2]} where 
   \eqn{n} is the number of parameters in par. Choices that exceed
   \eqn{2*n+1} are 
   not recommended.
   If not defined (i.e., the default \code{NA} is used), it will be set to
	\eqn{\min(n * 2, n+2)}{min(n * 2, n+2)}.
}
  \item{rhobeg}{\code{rhobeg} and \code{rhoend}
    must be set to the initial and final values of a trust
    region radius, so both must be positive with
    \code{0 < rhoend < rhobeg}. Typically
    \code{rhobeg} should be about one tenth of the greatest expected
    change to a variable.
}
\item{rhoend}{   The smallest value of the trust region radius that is allowed.
  If not defined, then 
   1e-6 times the value set for \code{rhobeg} will be used.
   \code{rhoend} should indicate the accuracy that is required in
   the final values of the variables.
}
  \item{iprint}{ The value of \code{iprint} should be set to \code{0, 1, 2 or 3},
  which controls the
  amount of printing. Specifically, there is no output if \code{iprint=0} and
       there is output only at the return if \code{iprint=1}. Otherwise, each
       new
       value of the trust region radius
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
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list;  see \code{\link{newuoa.control}} for
    details.
}

\keyword{optimize}

