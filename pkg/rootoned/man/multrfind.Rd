\name{multrfind}
\alias{multrfind}
\encoding{UTF-8}
\title{multrfind: Find a single root of a function of one variable within a specified interval.}
\concept{root-finding}
\description{
Find a single root of a function of one variable within a specified interval. Brent's method.}
\usage{
   multrfind(fn, gr=NULL, ri, ftrace=FALSE, meths="uniroot") 
}
\arguments{
 \item{fn}{A real function of a single variable.}
 \item{gr}{A real function of a single variable to evaluate the gradient of \code{fn}.}
 \item{ri}{A vector of two real numbers defining the ends of the interval in which to search.}
 \item{ftrace}{TRUE if the progress of the methods is to be displayed}
 \item{meths}{A vector of character (quoted) names of methods to apply to find real roots of \code{fn}.}
}
\details{
   This routine calls the rootfinder wrapper \code{rootwrap} for each of
   the named solvers. 
}

\value{
  The output is a data frame consisting of one row for each of the
  methods. The columns of the data frame are
  "root", "froot", "rtol", "iter", "fncount", and "method".
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

