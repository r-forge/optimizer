\name{newt1d}
\alias{newt1d}
\encoding{UTF-8}
\title{Find a single root of a function of one variable from an initial guess by Newton's method in 
    one dimension.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval.}
\usage{
    newt1d(fn, gr, x0, offset = 1000.0, maxiter = 100, trace=FALSE,...)
}
\arguments{
 \item{fn}{A real function of a single variable.}
 \item{gr}{A function that computes the gradient of \code{fn}.}
 \item{x0}{An initial guess for the value of the root of \code{fn}.}
 \item{offset}{An inverse tolerance for the root. When a new estimate of the root \code{xnew} is such 
    that \code{(xnew + offset)} is equal to \code{(xold+offset)}, \code{newt1d} will terminate.}
 \item{maxiter}{A limit on the number of function evaluations permitted in the attempt to find a root.}
 \item{trace}{If TRUE, provided diagnostic output of the process of rootfinding.}
 \item{...}{Extra data needed to compute the function \code{f}.}
}
\details{
   To be added.
}

\value{
  The output is a list consisting of 
    return(list(root=b, froot=NA, rtol=wtol, maxit=maxit)) ## Acceptable approx. is found	*/



  \item{root}{The estimate, b, of the position of the root.}
  \item{froot}{The value of the function at this point, f(b).}
  \item{rtol}{The width of the interval in which the root exists.}
  \item{fcount}{The number of function evaluations used.}
}
\seealso{
    \code{\link{optimize}}
}
\examples{
  # examples are in rootoned-package.Rd
}
\keyword{root-finding}

