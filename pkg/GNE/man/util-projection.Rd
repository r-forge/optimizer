\name{projector}
\alias{projector}
\alias{Projector}
\title{Projection of a point on a set}
\description{
		Projection of a point \code{z} on the set defined by the constraints \code{g(x) <= 0}.
}
\usage{
projector(z, g, jacg, bounds=c(0, 10), echo=FALSE, ...)
}
\arguments{
  \item{z}{The point to project.}
  \item{g}{The constraint function.}
  \item{jacg}{The jacobian of the constraint function.}
  \item{bounds}{bounds for the randomized initial iterate.}
  \item{echo}{a logical to plot traces.}
  \item{...}{further arguments to pass to \code{g} function.}
}
\details{
Find a point \code{x} in the set \eqn{K} which minimizes the Euclidean distance \eqn{||z - x||^2},
where the set \eqn{K} is \eqn{x, g(x) <= 0}. The Optimization is carried out by the \code{constrOptim.nl}
function of the package alabama.

}
\value{
A vector \code{x}. 
}
\seealso{
  See also \code{\link{GNE}}.
}
\author{
 Christophe Dutang
}
\examples{

# 1. the rectangle set
#

g <- function(x)
	c(x - 3, 1 - x)

jacg <- function(x)
	rbind(
	diag( rep(1, length(x)) ),
	diag( rep(-1, length(x)) )
	)

z <- runif(2, 3, 4)

#computation
projz <- projector(z, g, jacg)

#plot
plot(c(1, 3), c(1, 1), xlim=c(0, 4), ylim=c(0,4), type="l", col="blue")
lines(c(3, 3), c(1, 3), col="blue")
lines(c(3, 1), c(3, 3), col="blue")
lines(c(1, 1), c(3, 1), col="blue")

points(z[1], z[2], col="red")
points(projz[1], projz[2], col="red", pch="+")

z <- runif(2) + c(1, 0)
projz <- projector(z, g, jacg)

points(z[1], z[2], col="green")
points(projz[1], projz[2], col="green", pch="+")



# 2. the circle set
#

g <- function(x) sum((x-2)^2)-1
jacg <- function(x) as.matrix( 2*(x-2) )

z <- runif(2) + c(1, 0)

#computation
projz <- projector(z, g, jacg)

#plot
plot(c(1, 3), c(1, 1), xlim=c(0, 4), ylim=c(0,4), type="n", col="blue")
symbols(2, 2, circles=1, fg="blue", add=TRUE, inches=FALSE)

points(z[1], z[2], col="red")
points(projz[1], projz[2], col="red", pch="+")

z <- c(runif(1, 3, 4), runif(1, 1, 2))
projz <- projector(z, g, jacg)

points(z[1], z[2], col="green")
points(projz[1], projz[2], col="green", pch="+")


}
\keyword{nonlinear}
\keyword{optimize}

