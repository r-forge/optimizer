\name{polyarea}
\alias{polyarea}
\encoding{UTF-8}
\title{Area of a Polygon Given as Special Radial Coordinates.}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the area of the polygon.
}
\usage{
   polyarea(b)
}
\arguments{
 \item{b}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
}
\details{
   This routine computes area using the trigonometric formulas and does not
   first compute x, y coordinates of the vertices of the polygon. 
}
\value{
  The area of the polygon.
}

