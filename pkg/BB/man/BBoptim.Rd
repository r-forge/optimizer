\name{BBoptim}
\alias{BBoptim}
\title{Optimization Wrapper}
\description{A strategy to try different methods to obtain successful optimization.}
\usage{
  BBoptim(par, fn, gr=NULL, lower=-Inf, upper=Inf, method=c(2,3,1), 
        control=list(), ...) 
  }
\arguments{
    \item{par}{A real vector argument to \code{fn}, indicating the initial guess for the optimization of nonlinear 
    objective function \code{fn}.}

    \item{fn}{Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point.}
    \item{gr}{to add}
    \item{lower}{to add}
    \item{upper}{to add}
    \item{method}{A vector of integers specifying WHAT} 
    \item{control}{to add}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{NOT A list with elements as returned by \code{sane} or \code{dfsane}.}

\details{
Algorithm (\code{sane} or \code{dfsane}) is said to have failed when a 
non-zero convergence type is obtained, i.e. when \code{convergence > 0}. A user
might typically attempt various strategies, some of which are implemented in 
this wrapper. The function BBoptim tries the following sequential strategy:
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

\keyword{multivariate}
