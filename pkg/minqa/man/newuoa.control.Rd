\name{newuoa.control}
\Rdversion{1.1}
\alias{newuoa.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Controls for an R interface to the NEWUOA implementation of Powell
}
\description{
  The purpose of \code{newuoa.control} is to set control parameters for
  the NEWUOA Fortran code of M J D Powell. 
}
\usage{
newuoa.control(npt = NA, rhobeg = 0.01, rhoend = 1e-04, iprint = 3, maxfun = 10000, wsize = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npt}{
%%     ~~Describe \code{npt} here~~
   The number of points used to approximate the objective function via
   a quadratic approximation. Do we have a recommended setting?
}
  \item{rhobeg}{
%%     ~~Describe \code{rhobeg} here~~
}
  \item{rhoend}{
%%     ~~Describe \code{rhoend} here~~
}
  \item{iprint}{
%%     ~~Describe \code{iprint} here~~
  Output control -- integer -- with higher values yielding more output.
  Do we have a default that makes sense?
}
  \item{maxfun}{
%%     ~~Describe \code{maxfun} here~~
  The maximum allowed number of function evaluations. If this is exceeded, the
  method will terminate.
}  \item{wsize}{
%%     ~~Describe \code{wsize} here~~
  Allowed number of elements of working storage. Note that this should be dynamically
  allocated. It is NOT at present.
}
}
\details{
  Functions \code{fn} must return a numeric value.\cr
  The \code{control} argument is a list;  see \code{\link{newuoa.control}} for
    details.
}



\keyword{optimize}

