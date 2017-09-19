\name{addplot}
\alias{addplot}
\encoding{UTF-8}
\title{Save Polygon Parameters in Named Environment for Replotting.}
\description{
  Save current parameters, objective and area in named environment array \code{psave}
  and draw selected polygons.
}
\usage{
   addplot(penv, x, f, area) 
}
\arguments{
 \item{penv}{The (named) environment to be used to save polygon information.}
 \item{x}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
 \item{f}{The objective function computed for this set of parameters (given a particular
       \code{penfactor} and \code{epsilon}.}
 \item{area}{The polygon area for these parameters.}
}
\details{
   Draws the first (regular) polygon in pink and latest few in decreasing intensity
   (latest most dense). Saves point, area and function, and updates best area in
   \code{penv$besta}.
}
\value{
  0 for success, any other value for failure.
}

