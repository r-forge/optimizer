\name{hesschk}
\alias{hesschk}

\title{Run tests, where possible, on user objective function and (optionally) gradient and hessian}

\description{
   \code{hesschk} checks a user-provided R function, \code{ffn}. 
}

\usage{
   hesschk(xpar, ffn, ggr, hhess, trace=0, testtol=(.Machine$double.eps)^(1/3), ...)
}


\arguments{
     
    \item{xpar}{
        parameters to the user objective and gradient functions ffn and ggr
        }
    \item{ffn}{
        User-supplied objective function
        }
    \item{ggr}{
        User-supplied gradient function
        }
    \item{hhess}{
        User-supplied Hessian function
        }
    \item{trace}{
        set >0 to provide output from grchk to the console, 0 otherwise
        }
    \item{testtol}{
        tolerance for equality tests
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }

}


\details{
\tabular{ll}{
Package: \tab hesschk\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
\code{numDeriv} is used to compute a numerical approximation to the Hessian
matrix. If there is no analytic gradient, then the \code{hessian()} function 
from \code{numDeriv} is applied to the user function \code{ffn}. Otherwise, 
the \code{jacobian()} function of \code{numDeriv} is applied to the \code{ggr}
function so that only one level of differencing is used.
}


\value{
The function returns a single object \code{hessOK} which is TRUE if the 
analytic Hessian code returns a Hessian matrix that is "close" to the 
numerical approximation obtained via \code{numDeriv}; FALSE otherwise.

\code{hessOK} is returned with the following attributes:
\itemize{
     \item{"nullhess"}{Set TRUE if the user does not supply a function to compute the Hessian.}
     \item{"asym"}{Set TRUE if the Hessian does not satisfy symmetry conditions to within a tolerance.
              See the \code{hesschk} for details.}
     \item{"ha"}{The analytic Hessian computed at paramters \code{xpar} using \code{hhess}.}
     \item{"hn"}{The numerical approximation to the Hessian computed at paramters \code{xpar}.}
     \item{"msg"}{A text comment on the outcome of the tests.}
}


}


\author{

    John C. Nash
}

\examples{
# Want examples of both types of approximation to Hessian and also of failures
require("optfntools")

}

\keyword{optimize}

