\name{polysetup}
\alias{polysetup}
\encoding{UTF-8}
\title{Create Initial Polygon for R Optimization Animation}
\description{
  The purpose of \code{polysetup} is to create a polygon in the radial coordinate
  specification described in the Roptanimation vignette. For a polygon with nv
  vertices there are (nv - 1) radial parameters and (nv - 2) angular ones.
}
\usage{
   polysetup(nv, defsize=0.98)
}
\arguments{
 \item{nv}{The number of vertices in the polygon.}
 \item{defsize}{The default scale of the polygon to be created. Optimization methods
   generally do better with a polygon of scale smaller than 1, and the default for 
   \code{defsize} is 0.98.}
}
\details{
   None

}
\value{
  A list with components:
  \item{par0}{The set of radial and angular parameters describing a regular
    polygon of maximum diameter \code{defsize}.}
  \item{lb}{Lower bounds on these parameters -- all zeros.}
  \item{ub}{Upper bounds on these parameters, namely, 1 for the first (nv - 1) 
    parameters, and \code{pi} for the last (nv - 2). These are quite loose
    bounds for the angles.}
}
