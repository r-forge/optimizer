\name{uobyqa.control}
\Rdversion{1.1}
\alias{uobyqa.control}
\title{
Controls for an R interface to the uobyqa implementation of Powell
}
\description{
  The purpose of \code{uobyqa.control} is to set control parameters for
  the uobyqa Fortran code of M. J. D. Powell. 
}
\usage{
uobyqa.control(rhobeg = NA, rhoend = NA, iprint = 0, maxfun=10000)
}

\arguments{
  \item{rhobeg}{\code{rhobeg} and \code{rhoend}
    must be set to the initial and final values of a trust
    region radius, so both must be positive with
      \code{0 < rhoend < rhobeg}. Typically
    \code{rhobeg} should be about one tenth of the greatest expected
    change to a variable.}
  \item{rhoend}{
   The smallest value of the trust region radius that is allowed. If not defined, then 
   1e-6 times the value set for \code{rhobeg} will be used.
   \code{rhoend} should indicate the accuracy that is required in
   the final values of the variables.
}
\item{iprint}{
  The value of \code{iprint} should be set to \code{0, 1, 2 or 3},
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
  The maximum number of function evaluations allowed.
}
}
\details{
  Functions \code{fn} must return a numeric value.\cr
  The \code{control} argument is a list;  see \code{\link{uobyqa.control}} for
  details.  \cr

}




\keyword{optimize}

