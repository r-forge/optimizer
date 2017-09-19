\name{polygradu}
\alias{polygradu}
\encoding{UTF-8}
\title{Gradient of the log-barrier penalized negative area.}
\description{
  Using as input the particular radial coordinate specification given in the vignette, 
  compute the gradient of the log-barrier penalized negative area of the polygon.
}
\usage{
   polygradu(x, penfactor=1e-8, epsilon=0, penv)
}
\arguments{
 \item{x}{For a polygon in nv vertices, a vector of 2*nv - 3 parameters giving the radial
  and angular parameters that describe a particular polygon.}
 \item{penfactor}{Penalty factor for the log-barrier. Default 1E-8}
 \item{epsilon}{Shift to position of the barrier. Default 0}
 \item{penv}{Environment to save point information. (Unused here.)}
}
\details{
    The gradient of the penalized objective function is computed.   
}
\value{
  The vector that for the gradient is returned.
}