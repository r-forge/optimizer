\name{rejection}
\alias{rejection}
\alias{Rejection}
\title{Rejection method for random generation.}
\description{
		Generate random variate satisfying the constraint function by the Rejection algorithm.
}
\usage{
rejection(constr, nvars, LB=0, UB=1, ..., echo=FALSE, 
	method=c("unif","norm", "normcap"), control=list())

}
\arguments{
  \item{constr}{Constraint function}
  \item{nvars}{Number of variables}
  \item{LB}{Lower bound}
  \item{UB}{Upper bound}
  \item{...}{further arguments to pass to \code{constr} function.}
  \item{echo}{a logical to plot traces.}
  \item{method}{the distribution to draw random variates, either \code{"unif"},
	\code{"norm"}, \code{"normcap"}.}
}
\details{
Draw random variates \code{x} until all the components of \code{constr(x)} are negative. The distribution
to draw random variates can be the uniform distribution on the hypercube defined by \code{LB} and \code{UB},
the normal distribution centered in \code{(LB + UB)/2} and standard deviation \code{(UB - LB) / (4*1.9600)} 
and the capped normal distribution (intended for debug use).
}
\value{
A vector \code{x} which verifies the constraints  \code{constr(x) <= 0}.
}
\seealso{
  See also \code{\link{GNE}}.
}
\author{
 Christophe Dutang
}
\examples{

f <- function(x) x[1]^2 + x[2]^2 - 1

rejection(f, 2, -3, 3, method="unif")

rejection(f, 2, -3, 3, method="norm")

}
\keyword{optimize}

