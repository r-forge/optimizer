\name{polypar2XY}
\alias{polypar2XY}
\encoding{UTF-8}
\title{Compute x, y Coordinates for a Polygon Given in Particular Radial Coordinates}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the nv pairs of x,y coordinates, and append the first pair so that a single
  \code{plot} command gives a closed figure.
}
\usage{
   polypar2XY(b)
}
\arguments{
 \item{b}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
}
\details{
   None
}
\value{
  A list of 2 vectors. If the result is XY, then XY$x has (nv + 1) values for the x
  coordinates of the vertices (with the vertex 1 value repeated in first and last position),
  and XY$y has the corresponding y values.
}

