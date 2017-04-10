\name{snewton2}
\alias{snewton2}
\encoding{UTF-8}
\title{Safeguarded Newton method for function minimization using R functions.}
\description{This version of the safeguarded Newton solves the equations with
the R function solve() and uses the optimize() function as a line search.
}
\usage{
   snewton2(par, fn, gr, hess, control = list(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{gr}{A function that returns the gradient of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. This function 
   returns the gradient as a numeric vector.

 \item{hess}{A function to compute the Hessian matrix. This should be provided as a square,
      symmetric matrix.}
 \item{control}{
    An optional list of control settings.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value. \code{gr} must return a vector.
  \code{hess} must return a matrix. 
   The \code{control} argument is a list.
   \describe{
   \item{maxit}{A limit on the number of iterations (default 500 + 2*n where n is
      the number of parameters). This is the maximum number of gradient evaluations 
      allowed.}
   \item{maxfevals}{A limit on the number of function evaluations allowed 
     (default 3000 + 10*n).}
   \item{trace}{Set 0 (default) for no output, > 0 for diagnostic output
      (larger values imply more output).}
   \item{eps}{ a tolerance used for judging small gradient norm (default = 1e-07).
   	a gradient norm smaller than (1 + abs(fmin))*eps*eps is considered small 
   	enough that a local optimum has been found, where fmin is the current 
   	estimate of the minimal function value. }
   \item{acctol}{To adjust the acceptable point tolerance (default 0.0001) in the test
       ( f <= fmin + gradproj * steplength * acctol ). This test is used to ensure progress
       is made at each iteration. }
   \item{stepredn}{Step reduction factor for backtrack line search (default 0.2)}
   \item{reltest}{Additive shift for equality test (default 100.0)}
  }
}
}
\value{
  A list with components:
  \describe{
  \item{xs}{The best set of parameters found.}
  \item{fv}{The value of the objective at the best set of parameters found.}
  \item{grd}{The value of the gradient at the best set of parameters found. A vector.}
  \item{H}{The value of the Hessian at the best set of parameters found. A matrix.}
  \item{niter}{The number of Newton iterations used in finding the solution.}
  }
}
\references{ 

  Nash, J C (1979, 1990) Compact Numerical Methods for Computers: Linear
     Algebra and Function Minimisation, Bristol: Adam Hilger. Second
     Edition, Bristol: Institute of Physics Publications.

}
\seealso{\code{\link{optim}}}
\examples{
#Rosenbrock banana valley function
f <- function(x){
return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
#gradient
gr <- function(x){
return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
#Hessian
h <- function(x) {
a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
return(matrix(c(a11, a21, a21, 200), 2, 2))
}

fg <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  val
}
fgh <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  attr(val,"hessian") <- h(x)
  val
}

x0 <- c(-1.2, 1)

sol2 <- SNewton2(x0, fn=f, gr=gr, hess=h, control=list(trace=TRUE))
print(sol2)


}

\keyword{nonlinear}
\keyword{optimize}

