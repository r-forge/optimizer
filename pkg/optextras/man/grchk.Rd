\name{grchk}
\alias{grchk}

\title{Run tests, where possible, on user objective function and (optionally) gradient and hessian}

\description{
   \code{grchk} checks a user-provided R function, \code{ffn}. 
}


\usage{
   grchk(xpar, ffn, ggr, trace=0, testtol=(.Machine$double.eps)^(1/3), ...) 
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
Package: \tab grchk\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
\code{numDeriv} is used to numerically approximate the gradient of function \code{ffn}
and compare this to the result of function \code{ggr}.

}


\value{
\code{grchk} returns a single object \code{gradOK} which is true if the differences 
between analytic and approximated gradient are small as measured by the tolerance 
\code{testtol}.

This has attributes "ga" and "gn" for the analytic and numerically approximated gradients.

At the time of preparation, there are no checks for validity of the gradient code in
\code{ggr} as in the function \code{fnchk}.


}


\author{

    John C. Nash
}

\examples{
# Want examples of success and failure. What about "near misses"??


}

\keyword{optimize}

