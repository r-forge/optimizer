\name{internal}
\alias{squarem1}
\alias{squarem2}
\alias{cyclem1}
\alias{cyclem2}
\alias{poissmix.em}
\alias{poissmix.loglik}
\title{Internal functions for SQUAREM}
\description{SQUAREM functions not to be called by users.} 
\usage{
  poissmix.em(p,y)
  poissmix.loglik(p,y)
  squarem1(par, fixptfn, objfn, ... , control=list()) 
  squarem2(par, fixptfn, ... , control=list() )
  cyclem1(par, fixptfn, objfn, control=list(), ...)
  cyclem2(par, fixptfn, control=list(), ...)
  }
\arguments{
    \item{par}{See \code{\link{squarem}}.}
    \item{p}{See \code{\link{squarem}}.}
    \item{y}{See \code{\link{squarem}}.}
    \item{fixptfn}{See \code{\link{squarem}}.} 
    \item{objfn}{See \code{\link{squarem}}.} 
    \item{control}{See \code{\link{squarem}}.}
    \item{...}{See \code{\link{squarem}}.}
}
\keyword{internal}
