% File src/library/stats/man/zeroin.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{zeroin}
\alias{zeroin}
\encoding{UTF-8}
\title{zeroin: Find a single root of a function of one variable within a specified interval.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval.}
\usage{
    zeroin(f, interval,tol = .Machine$double.eps^0.5, maxiter = 1000, trace=FALSE, ...) 
}
\arguments{
 \item{f}{A real function of a single variable.}
 \item{interval}{A vector of two real numbers defining the ends of the interval in which to search.}
 \item{tol}{A tolerance for the root. When the interval within which a root exists is smaller than
 \code{tol}, zeroin will terminate.}
 \item{maxiter}{A limit on the number of function evaluations permitted 
 in the attempt to find a root.}
 \item{trace}{If TRUE, provided diagnostic output of the process of rootfinding.}
 \item{...}{Extra data needed to compute the function \code{f}.}
}
\details{
   This routine was created at UseR!2001 at the University of Warwick during a session 
   where there were no presentations of particular interest to the author (JN). It
   is based on \code{zeroin.c} (from the R source code) and \code{zeroin.f} (from
   the Netlib Forsythe, Malcolm and Moler collection at  https://netlib.sandia.gov/fmm/).
   
   This routine is entirely in \code{R} so that it can be used within package \code{Rmpfr}
   which allows long-precision numbers and arithmetic. Martin Maechler, author of 
   \code{Rmpfr}, adapted the present code to his package where it is called \code{unirootR}.
}


\value{
  The output is a list consisting of 
\itemize{ 
  \item root The estimate, b, of the position of the root.
  \item froot The value of the function at this point, f(b).
  \item rtol  The width of the interval in which the root exists.
  \item fcount  The number of function evaluations used.
}
}
\references{
  Brent, R. (1973)
  \emph{Algorithms for Minimization without Derivatives.}
  Englewood Cliffs, NJ: Prentice-Hall.
}

\seealso{
    \code{\link{optimize}}
}
\examples{
# Dekker example
require(rootoned)
dek <- function(x){ 1/(x-3) - 6 }
r1 <- zeroin(dek, c(3.0000001, 6), trace=TRUE)
r1
}

\keyword{root-finding}

