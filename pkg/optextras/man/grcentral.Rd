\name{grcentral}
\alias{grcentral}

\title{Central difference numerical gradient approximation.}

\description{
   \code{grcentral} computes the central difference approximation to the gradient of 
          user function \code{userfn}. 
}


\usage{
   grcentral(par, userfn, eps=1.0e-7, ...)
}


\arguments{
     
    \item{par}{
        parameters to the user objective function userfn
        }
    \item{userfn}{
        User-supplied objective function
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
Package: \tab grcentral\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
}


\value{
\code{grcentral} returns a single vector object \code{df} which approximates the 
gradient of userfn at the parameters par.

}


\author{

    John C. Nash
}

\examples{
cat("Example of use of grfwd\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grcentral(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)

}

\keyword{optimize}

