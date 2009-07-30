\name{bobyqa.control}
\Rdversion{1.1}
\alias{bobyqa.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Controls for an R interface to the bobyqa implementation of Powell
}
\description{
  The purpose of \code{bobyqa.control} is to set control parameters for
  the bobyqa Fortran code of M J D Powell. 
}
\usage{
%- bobyqa.control(npt = NA, rhobeg = 0.01, rhoend = 1e-04, iprint = 3, maxfun = 10000, wsize = 10000)
bobyqa.control(npt = NA, rhobeg = NA, rhoend = NA, iprint = 0,
maxfun = 10000, wsize = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npt}{
%%     ~~Describe \code{npt} here~~
   The number of points used to approximate the objective function via
   a quadratic approximation. The value of npt must be in the interval 
   [n+2,(n+1)(n+2)/2] where 
   n is the number of parameters in par. Choices that exceed 2*n+1 are 
   not recommended.
   If not defined (i.e., the default NA is used), it will be set to
	min(n * 2, n+2)
}
  \item{rhobeg}{
%%     ~~Describe \code{rhobeg} here~~
   The quadratic approximations estimate the gradient of the function
   at the current estimate of the minimum. The approximation needs the
   estimating points to be "not too close" together. A lower bound to
   this measure of distance is rho, and this is initially set to 
   \code{rhobeg}. If the user does not provide a value, this will be set
   to the maximum box constraint interval (for finite limits), i.e.,
        max(xu-xl)
   or to 
        abs(max(par) / 2)
   if the bounds are infinite (i.e., no box constraints are specified)
}
  \item{rhoend}{
%%     ~~Describe \code{rhoend} here~~
   The smallest value of rho that is allowed. If not defined, then 
   1e-6 times the value set for \code{rhobeg} will be used.
   \code{rhobeg} and \code{rhoend} must be positive with rhoend no 
   greater than rhobeg. Otherwise bobyqa will stop.
}
  \item{iprint}{
%%     ~~Describe \code{iprint} here~~
  Output control -- integer -- with higher values yielding more output.
  Default value is 0 (no output).
}
  \item{maxfun}{
%%     ~~Describe \code{maxfun} here~~
  The maximum allowed number of function evaluations. If this is exceeded, the
  method will terminate.
}  \item{wsize}{
%%     ~~Describe \code{wsize} here~~
  Allowed number of elements of working storage. If this is undefined (i.e., 
  left at the default NA), then it will be set to 
      (npt+5)*(npt+n)+3*n*(n+5)/2
  Note that npt needs to be set (and will be by default) BEFORE wsize can be
  set to this minimum size.
}
}
\details{
  Functions \code{fn} must return a numeric value.\cr
  The \code{control} argument is a list;  see \code{\link{bobyqa.control}} for
    details.
}

\keyword{optimize}

