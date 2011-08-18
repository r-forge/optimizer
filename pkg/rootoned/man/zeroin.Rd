% File src/library/stats/man/zeroin.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{zeroin}
\alias{zeroin}
\encoding{UTF-8}
\title{Find a single root of a function of one variable within a specified interval.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval.}
\usage{
    zeroin(f, interval,tol = .Machine$double.eps^0.5, maxiter = 1000, ...) 
}
\arguments{
 \item{f}{A real function of a single variable.}
 \item{interval}{A vector of two real numbers defining the ends of the interval in which to search.}
 \item{tol}{A tolerance for the root. When the interval within which a root exists is smaller than \code{tol}, 
    zeroin will terminate.}
 \item{maxiter}{A limit on the number of function evaluations permitted in the attempt to find a root.}
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

