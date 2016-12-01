\name{polypardist2}
\alias{polypardist2}
\encoding{UTF-8}
\title{Vertex to Vertex Distances for a Polygon Given as Special Radial Coordinates.}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the (nv - 1) * nv / 2 distances between all the vertices.
}
\usage{
   polypardist2(b)
}
\arguments{
 \item{b}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
}
\details{
   This routine computes distances using the trigonometric formulas and does not
   first compute x, y coordinates of the vertices of the polygon. 
}
\value{
  A vector of the distances in order (1-2), (1-3), ..., (1-(nv-1)), (2-3), ...
}
\seealso{\code{polypar2distXY}}

