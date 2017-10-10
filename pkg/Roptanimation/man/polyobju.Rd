\name{polyobju}
\alias{polyobju}
\encoding{UTF-8}
\title{Area of a Polygon Given as Special Radial Coordinates.}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the objective function that is negative area plus penalized sum of the 
  logs of slacks in the vertex distance constraints adjusted by epsilon.
  This version provides barrier constraints for the radial parameters as well as
  non-radial distances.
}
\usage{
   polyobju(x, penfactor=1e-5, epsilon=0, penv)
}
\arguments{
 \item{x}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
 \item{penfactor}{Penalty factor for the log-barrier. Default 1E-5}
 \item{epsilon}{Shift to position of the barrier. Default 0}
 \item{penv}{Environment where point data is to be saved.}
}
\details{
    The negative area is computed via function polyarea(). 

    The slacks are, for each vertex-vertex distance, given by
       slack <- 1 + epsilon - distance

    The negative area + penfactor * sum(log(slacks)) is the objective function.
   
}
\value{
  The log-barrier penalized objective function with attributes "area" and 
  "minslack" = min(slacks).
}

