\name{rootwrap}
\alias{rootwrap}
\encoding{UTF-8}
\title{zeroin: Find a single root of a function of one variable within a specified interval.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval. Brent's method.}
\usage{
   rootwrap(fn=fn, gr=NULL, ri=NULL, method="uniroot", ftrace=TRUE, ...)
}
\arguments{
 \item{fn}{A real function of a single variable for which a root is desired.}
 \item{gr}{A real function of a single variable giving the gradient (derivative) of \code{fn}.}
 \item{ri}{A vector of two real numbers defining the ends of the interval in which to search.}
 \item{method}{A character string giving the name of the rootfinding method to be used. Must be
       one of the methods \code{rootwrap} supports. See the source code.}
 \item{ftrace}{If TRUE, provided diagnostic output of the process of rootfinding.}
 \item{...}{Extra data needed to compute the function \code{f}.}
}
\details{
    This routine is a wrapper to allow a consistent calling sequence for a number of different
    rootfiners. 
}

\value{
  The output is a list consisting of 
\itemize{  
  \item{root} {The estimate, b, of the position of the root.}
  \item{froot} {The value of the function at this point, f(b).}
  \item{rtol} {The width of the interval in which the root exists.}
  \item{iter} {The number of iteration cycles used.}
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
# require(rootoned)
dek <- function(x){ 1/(x-3) - 6 }
r1 <- rootwrap(dek, ri=c(3.0000001, 6), ftrace=TRUE, method="uniroot")
r1
}
\keyword{root-finding}

