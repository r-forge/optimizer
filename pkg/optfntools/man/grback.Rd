\name{grback}
\alias{grback}

\title{Backward difference numerical gradient approximation.}

\description{
   \code{grback} computes the backward difference approximation to the gradient of 
          user function \code{userfn}. 
}


\usage{
   grback(par, userfn, fbase=NULL, eps=1.0e-7, ...)
}


\arguments{
     
    \item{par}{
        parameters to the user objective function userfn
        }
    \item{userfn}{
        User-supplied objective function
        }
    \item{fbase}{
        The value of the function at the parameters, else NULL. This is to save
        recomputing the function at this point.
        }
    \item{eps}{
        a small number used to generate the finite difference step
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}


\details{
\tabular{ll}{
Package: \tab grback\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
}


\value{
\code{grback} returns a single vector object \code{df} which approximates the 
gradient of userfn at the parameters par.

}


\author{

    John C. Nash
}

\examples{
# Want examples and comparisons with numDeriv and other forms. ??

require("numDeriv")





}

\keyword{optimize}

