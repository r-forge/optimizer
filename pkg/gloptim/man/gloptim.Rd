\name{gloptim}
\alias{gloptim}
\title{
  Wrapper to unify global optimization solvers
}
\description{
  A routine to call a variety of global optimization methods through a
  single syntax. The structure is intended to be extensible.
}
\usage{
gloptim(fn, lb, ub, x0 = NULL,
        method = c("deoptim", "deoptimr", "simplede", "ga", "smco", "soma"),
        type = NULL,
        minimize = TRUE, control = list(), ...)
}

\arguments{
  \item{fn}{Nonlinear objective function that is to be optimized. 
            A scalar function that takes a real vector as argument and 
            returns a scalar that is the value of the function at that point.}
  \item{lb}{Lower bounds on the parameters.
            Vector the same length as the variables of the objective function.}
  \item{ub}{Upper bounds on the parameters.
            Vector the same length as the variables of the objective function.}
  \item{x0}{Starting vector of parameter values.
            For some methods, this is optional.}
  \item{method}{A character string specifying the solver to be applied.
                Current options are "deoptim[r]", "ga", "simplede", "soma", and "smco".}
  \item{type}{A characterization of the function that may be needed by some methods.}
  \item{minimize}{TRUE (default) if the function is to be minimized.}
  \item{control}{A list of control parameters.
       See \bold{Details} for more information.}
  \item{\dots}{Additional arguments passed to the objective function \code{fn}.}
}
\details{
  Argument \code{control} is a list specifing changes to default 
  values of algorithm control parameters.

  If the minimization process threatens to go into an infinite loop,
  try to set \code{maxiter}.
}
\value{
  A list with the following components:
  \item{xmin}{Best estimate of the parameter vector found by the algorithm.}
  \item{fmin}{value of the objective function at termination.}
}
\author{
  Hans W Borchers  <hwborchers@googlemail.com> and \cr
  John C. Nash <nashjc@uottawa.ca>
}
\seealso{
  \code{\link{optim}}
}
\examples{
##  See the tests directory and the vignettes.
}
\keyword{ global optimization }
\keyword{ stochastic optimization }

