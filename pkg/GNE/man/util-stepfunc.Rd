\name{stepfunc}
\alias{stepfunc}
\alias{purestep}
\alias{decrstep}
\alias{decrstep10}
\alias{decrstep20}
\alias{decrstep30}
\title{Step functions}
\description{
		Step functions for relaxation methods
	}
\usage{
purestep(k)
decrstep(k, param)
decrstep10(k)
decrstep20(k)
decrstep30(k)
}
\arguments{
  \item{k}{iteration number.}
  \item{param}{parameter for the decreasing step function after which the step decreases.}
}
\details{
The \code{decrstep} function is a decreasing step serie such that \code{decrstep(k)} 
equals to \eqn{1/(k - param)}. Functions \code{decrstep10, decrstep20, decrstep30} are just wrappers of
\code{decrstep}. 

The \code{purestep} function implements a constant step serie equaled to 1.
}
\value{
A numeric.
}
\seealso{
  See also \code{\link{GNE}}.
}
\author{
 Christophe Dutang
}
\examples{

purestep(13)

decrstep(13, 7)

decrstep10(13)

decrstep20(13)

decrstep30(13)

}
\keyword{math}


