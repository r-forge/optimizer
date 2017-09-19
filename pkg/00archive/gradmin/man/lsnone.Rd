\name{lsnone}
\alias{lsnone}
\encoding{UTF-8}
\title{Null line search for gradmin.}
\description{A "line search" that always takes a unit step.
}
\usage{
   lsnone<-function(fn, fbest, xc, d, grv, ws, \dots)
}
\arguments{
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{fbest}{Current lowest function value found.}
 \item{xc}{parameter vector that gives lowest function value so far}
 \item{d}{search direction from point \code{xc}.}
 \item{grv}{gradient vector at the point \code{xc}}
 \item{ws}{A workspace, in the form of a list that is passed to other functions. This 
      workspace can be used to hold and pass various computed elements or control parameters.}
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}

\details{
  Functions \code{fn} must return a numeric value. \code{gr} must return a vector.
}

\value{
  \code{st}, a step from \code{xc} along \code{d} that satisfies the Armijo
  condition 
  
    \code{f(xc + st*d) <= fbest + acctol*st*gproj}
    
    where \code{acctol} is a tolerance (provided in the \code{control} vector)
    and \code{grproj} is the projection of the gradient \code{grv} along
    \code{d}. (This projection should be negative.)

  The attribute \code{Fval} of the step is the function value for that step
  size, that is, at \code{xc + st*d}
}
\references{ 

  Nash, J C (1979, 1990) Compact Numerical Methods for Computers: Linear
     Algebra and Function Minimisation, Bristol: Adam Hilger. Second
     Edition, Bristol: Institute of Physics Publications.

}
\seealso{\code{\link{optim}}}
\keyword{nonlinear}
\keyword{optimize}

