\name{stepfunc}
\alias{stepfunc}
\alias{purestep}
\alias{decrstep}
\alias{decrstep5}
\alias{decrstep10}
\alias{decrstep20}
\title{Step functions}
\description{
		Step functions for relaxation methods
	}
\usage{
purestep(k)
decrstep(k, param)
decrstep5(k)
decrstep10(k)
decrstep20(k)
}
\arguments{
  \item{k}{iteration number.}
  \item{param}{parameter for the decreasing step function after which the step decreases.}
}
\details{
The \code{decrstep} function is a decreasing step serie such that \code{decrstep(k)} 
equals to \eqn{1/2/(k - param)} when \eqn{k>param}, 1/2, otherwise. 
Functions \code{decrstep5, decrstep10, decrstep20} are just wrappers of
\code{decrstep}. 

The \code{purestep} function implements a constant step serie equaled to 1.
}
\value{
A numeric.
}
\seealso{
  See also \code{\link{GNE}} and \code{\link{GNE.fpeq}}.
}
\author{
 Christophe Dutang
}
\examples{

cbind(
purestep(1:20),
decrstep(1:20, 7),
decrstep5(1:20),
decrstep10(1:20),
decrstep20(1:20)
)

}
\keyword{math}


