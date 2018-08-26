\name{FnTrace}
\alias{FnTrace}
\alias{TraceSetup}
\alias{grTrace}
\encoding{UTF-8}
\title{FnTrace: A function wrapper for rootfinding to allow diagnostic tracking of progress.}
\concept{root-finding}
\description{
Track progress of rootfinding.}
\usage{

   TraceSetup(ifn=0, igr=0, ftrace=FALSE, fn=NA, gr=NA)
   
   FnTrace(x,...)
   
   grTrace(x,...)
}
\arguments{
 \item{ifn}{Count of the number of function evaluations.}
 \item{igr}{Count of the number of gradient evaluations.}
 \item{ftrace}{If TRUE, provide diagnostic output of the process of rootfinding.}
 \item{fn}{The function for which roots are to be found.}
 \item{gr}{The function giving the gradient of \code{gr}}
 \item{x}{If used, the argument for which \code{fn} and \code{gr} are to be evaluated.}
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
r1 <- root1d(dek, c(3.0000001, 6), trace=TRUE)
r1
}
\keyword{root-finding}

