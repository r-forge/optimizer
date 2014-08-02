\name{lbfgsb3}
\alias{lbfgsb3}
\title{Interfacing to Nocedal LBFGSB3 limited memory BFGS solver.
}
\description{
   Interfacing to Nocedal LBFGSB3 limited memory BFGS solver.
}
\usage{
   lbfgsb3(prm, fn, gr=NULL, lower=-Inf, upper=Inf, 
         control=list(), ...)
}
\arguments{
  \item{prm}{
     A parameter vector which gives the initial guesses to the parameters
     that will minimize \code{fn}. This can be named, for example, we could use
        prm=c(b1=1, b2=2.345, b3=0.123)
  }
  \item{fn}{
     A function that evaluates the objective function to be minimized.
  }
  \item{gr}{
     If present, a function that evaluates the gradient vector for 
     the objective function at the given parameterscomputing the elements of
     the sum of squares function at the set of parameters \code{start}. 
  }
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default -Inf.
  } 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default Inf.
  } 
  \item{control}{ 
     A list of controls for the algorithm. These are:
     \describe{
      \item{\code{femax}}{Maximum function evaluations. }
      \item{\code{gemax}}{Maximum number of gradient evaluations. }
     }
  }
  \item{\dots}{
     Any data needed for computation of the objective function and gradient.
}

}
\details{
    To be added.
}
\value{
  A list of the following items
  \item{prm}{A vector giving the parameter values at the supposed solution.}
  \item{f}{The value of the objective function at this set of parameters.}
  \item{g}{An estimate of the gradient of the objective at the solution.}
}
\references{
%% ~put references to the literature/web site here ~
     Nocedal et al. (1994) 

     others!!

}
\author{
John C Nash <nashjc@uottawa.ca>
}
\note{
   Special notes, if any, will appear here.
}

\seealso{
    Packages \code{\link{optim}} and \code{optimx}.
}

\examples{
cat("Examples are to be added\n")

}
\keyword{ nonlinear parameter optimization }

