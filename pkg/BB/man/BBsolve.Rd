\name{BBsolve}
\alias{BBsolve}
\title{Optimization Wrapper}
\description{A strategy to try different methods to obtain successful optimization.}
\usage{
  BBsolve(par, fn, algorithm = c("dfsane", "sane"), 
        method=c(2,3,1), M = c(10, 50), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500) 
  }
\arguments{
    \item{par}{A real vector argument to \code{fn}, indicating the initial guess for the optimization of nonlinear 
    objective function \code{fn}.}

    \item{fn}{Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point.}
    \item{algorithm}{character string indicating "sane" or "dfsane".} 
    \item{method}{A vector of integers specifying WHAT} 
    \item{M}{to add} 
    \item{NM}{to add}
    \item{tol}{to add}
    \item{maxit}{to add}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{NOT A list with elements as returned by \code{sane} or \code{dfsane}.}

\details{
Algorithm (\code{sane} or \code{dfsane}) is said to have failed when a 
non-zero convergence type is obtained, i.e. when \code{convergence > 0}. A user
might typically attempt various strategies, some of which are implemented in 
this wrapper. The function BBsolve tries the following sequential strategy:
\enumerate{
\item Try a different starting value. This should be tried first when sensible 
   starting values can be generated easily (e.g. randomly perturbing the 
   original starting value). 
\item Try a different method.  Since the default is \code{method = 2}, the 
   user should try methods 1 and/or 3. 
\item Try a different non-monotonicity parameter \code{M}.  Since the default 
   is \code{M = 10}, the user may try \code{M=100}. 
\item Try with Nelder-Mead initialization \code{NM}.  Since the default 
   is \code{NM = FALSE}, the user should try \code{NM = TRUE}. 
}
}
\seealso{
  \code{\link{sane}},
  \code{\link{dfsane}}
   }
   
\examples{
# Use a preset seed so test values are reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

broydt <- function(x) {
n <- length(x)
f <- rep(NA, n)
h <- 2
f[1] <- ((3 - h*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - h*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - h*x[n]) * x[n]) - x[n-1] + 1
f
}

BBsolve(par=rnorm(500), fn=broydt)  # note that the default doesn't work generally.
}

\keyword{multivariate}
