% File src/library/stats/man/root1d.Rd

\name{root1d}
\alias{root1d}
\encoding{UTF-8}
\title{root1d: Find a single root of a function of one variable within a specified interval.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval.}
\usage{
   root1d(f, interval, tol = .Machine$double.eps^0.5, maxiter = 1000, trace=FALSE, ...) 
}
\arguments{
 \item{f}{A real function of a single variable.}
 \item{interval}{A vector of two real numbers defining the ends of the interval in which to search.}
 \item{tol}{A tolerance for the root. When the interval within which a root exists is smaller than \code{tol}, 
    root1d will terminate.}
 \item{maxiter}{A limit on the number of function evaluations permitted in the attempt to find a root.}
 \item{trace}{If TRUE, provided diagnostic output of the process of rootfinding.}
 \item{...}{Extra data needed to compute the function \code{f}.}
}
\details{
   This routine is an all-R translation of Algorithm 18 of Nash (1979). This is a
   variant of the Dekker-Brent approach. At the time it was written, the goal was
   to produce a very short code that was robust in a low-precision environment. 
}

\value{
  The output is a list consisting of 
\itemize{  
  \item{root} {The estimate, b, of the position of the root.}
  \item{froot} {The value of the function at this point, f(b).}
  \item{rtol} {The width of the interval in which the root exists.}
  \item{fcount} {The number of function evaluations used.}
}  
}

\references{
  Nash, John C. (1979)
  \emph{Compact numerical methods: linear algebra and function minimisation.}
  Bristol UK: Adam Hilger.
}

\seealso{
    \code{\link{optimize}}
}
\examples{
# Dekker example
require(rootoned)
dek <- function(x){ 1/(x-3) - 6 }
r1 <- root1d(dek, c(3.0000001, 6), trace=TRUE)
r1
}
\keyword{root-finding}

