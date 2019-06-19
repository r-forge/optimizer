\name{checksolver}
\alias{checksolver}
\encoding{UTF-8}
\title{Test if requested solver is present}
\concept{minimization}
\concept{maximization}
\description{
  Test if requested solver is present.
}
\usage{
   checksolver(method, allmeth, allpkg)
}
\arguments{
 \item{method}{Character string giving the name of the solver requested.}
 \item{allmeth}{Character vector giving the names of the methods optimr can use.}
 \item{allpkg}{Character vector giving the names of the packages where
     solvers are found.}
}
\value{
  \code{checksolver} tests if requested function minimization solver is present.
}
\keyword{nonlinear}
\keyword{optimize}

