\name{polypar2distXY}
\alias{polypar2distXY}
\encoding{UTF-8}
\title{Vertex to Vertex Distances for a Polygon Given as Special Radial Coordinates.}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the (nv - 1) * nv / 2 distances between all the vertices.
}
\usage{
   polypar2distXY(pars)
}
\arguments{
 \item{pars}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
}
\details{
   This routine calls \code{polypar2XY} to get the x, y coordinates of the polygon,
   then calls \code{polydistXY} to get the distances.
}
\value{
  A vector of the distances in order (1-2), (1-3), ..., (1-(nv-1)), (2-3), ...
}
\seealso{\code{polypardist2}}