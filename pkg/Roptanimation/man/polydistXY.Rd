\name{polydistXY}
\alias{polydistXY}
\encoding{UTF-8}
\title{Vertex to Vertex Distances for a Polygon Given as XY Plotting Coordinates.}
\description{
  Using as input the x and y coordinates of nv + 1 points, where the last point is 
  the same as the first to allow for easy drawing of the closed polygon in a single
  \code{plot} command, compute the (nv - 1) * nv / 2 distances between all the 
  vertices.
}
\usage{
   polydistXY <- function(XY)
}
\arguments{
 \item{XY}{A list of two vectors, \code{XY$x} and \code{XY$y} giving the vertices of
  the polygon in order, with a repetition of the first x or y appended to each of the
  vectors.}
}
\details{
   None

}
\value{
  A vector of the distances in order (1-2), (1-3), ..., (1-(nv-1)), (2-3), ...
}
