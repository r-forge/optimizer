\name{multiStart}
\alias{multiStart}
\title{Optimization Starting Point Sensitivity Test}
\description{Start BBsolve from multiple starting points to check sensitivity.}
\usage{
  multiStart(par, fn, algorithm = c("dfsane", "sane"), 
	method=c(2,3,1), M = c(10, 50), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500, details=FALSE) 
  }
\arguments{
    \item{par}{A real matrix, each row of which is an argument to \code{fn}, 
       indicating initial guesses for the optimization of the objective 
       function \code{fn}.}
    \item{fn}{see BBsolve.} 
    \item{algorithm}{see BBsolve.} 
    \item{method}{see BBsolve.} 
    \item{M}{see BBsolve.} 
    \item{NM}{see BBsolve.}
    \item{tol}{see BBsolve.}
    \item{maxit}{see BBsolve.}
    \item{details}{logical indicating if \code{info} should be returned.}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{list with elements \code{par}, \code{values}, \code{converged}, 
  and optionally \code{info}. Each element of
  \code{info} is a list as returned by \code{BBsolve}, which in turn is a list
  as returned by \code{sane}or \code{dfsane}. The elements  \code{par}, 
  \code{values}, and \code{converged} are all information from \code{info} in
  a more convenient vector or matrix format}

\details{
The optimization is run with each row of \code{par} indicating initial guesses for 
for the optimization. If \code{details} is \code{TRUE} then the returned list 
includes an element \code{info} which is a list with each element the complete
returned value from the optimization for the corresponding \code{par} initial guess.
}
\seealso{
  \code{\link{BBsolve}},
  \code{\link{sane}},
  \code{\link{dfsane}}
   }
   
\examples{
# Use a preset seed so the example is reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

# Brown's almost linear system(A.P. Morgan, ACM 1983)
brownlin <- function(x) {
# two distinct solutions if n is even
# three distinct solutions if n is odd
  	n <- length(x)
  	f <- rep(NA, n)
	nm1 <- 1:(n-1)
	f[nm1] <- x[nm1] + sum(x) - (n+1)
	f[n] <- prod(x) - 1 
	f
}

p0 <- matrix(rnorm(200), 20, 10)  # 20 starting values, each of length 10
ans <- multiStart(par=p0, fn=brownlin)

pc <- princomp(ans$par) 
plot(pc$scores[,1])  # plot shows two/three distinct solutions depending on n is even/odd
}

\keyword{multivariate}
