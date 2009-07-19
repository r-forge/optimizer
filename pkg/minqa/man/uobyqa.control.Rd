\name{uobyqa.control}
\Rdversion{1.1}
\alias{uobyqa.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Controls for an R interface to the uobyqa implementation of Powell
}
\description{
  The purpose of \code{uobyqa.control} is to set control parameters for
  the uobyqa Fortran code of M J D Powell. 
}
\usage{
uobyqa.control(trace = 0, rhobeg =  log(100)/10, rhoend = 1e-04, maxit = 1e+05, 
     fnscale = 1, parscale = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trace}{
%%     ~~Describe \code{npt} here~~
   Output control. Larger numbers give more output. Default is 0. (No output.)
}
  \item{rhobeg}{
%%     ~~Describe \code{rhobeg} here~~
}
  \item{rhoend}{
%%     ~~Describe \code{rhoend} here~~
}
  \item{maxit}{
%%     ~~Describe \code{iprint} here~~
  The maximum number of 'iterations'. Do we know what this means?
}
  \item{fnscale}{
%%     ~~Describe \code{maxfun} here~~
  A multiplicative scaling factor for the objective function. Set this to
  -1 to maximize a function.}  
  \item{parscale}{
%%     ~~Describe \code{parscale} here~~
  A (vector) parameter scaling to try to render the scale of the parameters
  within the minimization algorithm approximately similar.}
  \item{\dots}{
    Any other controls that might be needed. Why is this here?
  }
}
\details{
  Functions \code{fn} must return a numeric value.\cr
  The \code{control} argument is a list;  see \code{\link{uobyqa.control}} for
    details.
}



\keyword{optimize}

