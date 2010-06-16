\name{squarem}
\alias{squarem}
\alias{squarem1}
\alias{squarem2}
\alias{poissmix}
\alias{poissmix.em}
\alias{poissmix.loglik}
\alias{cyclem1}
\alias{cyclem2}
\alias{emiter}
\title{(Squared extrapolation methods for optimization acceleration FIX}
\description{FIX}
\usage{
  poissmix.loglik(p,y)
  poissmix.em(p,y)
  poissmix.loglik(p,y)
  squarem(par, fixptfn, objfn, ... , control=list())
  squarem1(par, fixptfn, objfn, ... , control=list())
  squarem2(par, fixptfn, ... , control=list() )
  cyclem1(par, fixptfn, objfn, control=list(), ...)
  cyclem2(par, fixptfn, control=list(), ...)
  emiter(p, fixptfn, objfn, tol=1.e-07, maxiter=5000, trace=FALSE, ...)
  }
\arguments{
    \item{p}{FIX.} 
    \item{y}{FIX.} 
    \item{par}{FIX A real matrix, each row of which is an argument to \code{fn},
      indicating initial guesses for solving a nonlinear system \code{fn = 0} 
      or for optimizing the objective function \code{fn}.}
    \item{fixptfn}{FIX.} 
    \item{objfn}{FIX.} 
    \item{control}{FIX See \code{BBsolve} and \code{BBoptim}.}
    \item{tol}{FIX.}
    \item{maxiter}{FIX.}
    \item{trace}{FIX.}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{whatever.}

\details{
The optimization ... .
}
\seealso{
  \code{\link[BB]{BBsolve}},
  \code{\link[BB]{BBoptim}},
  \code{\link[BB]{dfsane}},
  \code{\link[BB]{spg}}
   }
   
\examples{
# Use a preset seed so the example is reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

# FIX add more


}
\references{ 
   FIX R Varadhan Scand J Stats EXAMPLEand PD Gilbert (2009),  BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, \url{http://www.jstatsoft.org/v32/i04/}
}

\keyword{multivariate}
