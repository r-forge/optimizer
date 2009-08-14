\name{Rcgmin}
\alias{Rcgmin}
\alias{Rcgminu}
\encoding{UTF-8}
\title{An R implementation of a nonlinear conjugate gradient algorithm with the Dai / Yuan
   update and restart. Based on Nash (1979) Algorithm 22 for main structure.}
\description{
  The purpose of \code{Rcgminu} is to minimize an unconstrained function 
  of many parameters by a nonlinear conjugate gradients method. The variant
  \code{Rcgmin} permits bounds and mask constraints to be imposed on the 
  parameters.
}
\usage{
   Rcgminu(par, fn, gr, control = list(), \dots)
   Rcgmin(par, fn, gr, lower, upper, bdmsk, control = list(), \dots)
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

   If gr is NOT supplied, a numerical approximation will be used. This 
   is STRONGLY discouraged, as the CG methods seem to require accurate
   gradients. Furthermore, the computation of numerical gradients is
   many times slower -- of the order of 1000s -- than using an analytic
   function.
}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization.}
 \item{control}{
    An optional list of control settings.  
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Functions \code{fn} must return a numeric value.
  The \code{control} argument is a list.
    Successful completion.
    The codes \code{Rcgmin} for R are still a work in progress, so users should watch
    the console output.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
}
\references{
To be added.
}
\seealso{\code{\link{optim}}}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
x3 <- Rcgminu(fn=fr,par=c(1,2))
}
\keyword{nonlinear}
\keyword{optimize}

