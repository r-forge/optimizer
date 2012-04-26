\name{modss}
\alias{modss}
\title{Compute gradient from residuals and Jacobian.
}
\description{
   For a nonlinear model originally expressed as an expression of the form
         lhs ~ formula_for_rhs
   assume we have a resfn and jacfn that compute the residuals and the 
   Jacobian at a set of parameters. This routine computes the gradient, 
   that is, t(Jacobian) %*% residuals. 
}
\usage{
   modss(prm, resfn, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prm}{
     A parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     However, the names are NOT used, only positions in the vector.
}
  \item{resfn}{
     A function to compute the residuals of our model at a parameter vector.
} 
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{modss} calls resfn to compute residuals and then uses \code{crossprod}
    to compute the sum of squares. 

    At 2012-4-26 there is no checking for errors. 

    Note that it appears awkward to use this function in calls to optimization
    routines. The author would like to learn why. 
}
\value{
    The numeric value of the sum of squares at the paramters.
}
\references{
%% ~put references to the literature/web site here ~
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!

}
\author{
John C Nash <nashjc@uottawa.ca>
}
\note{
   Special notes, if any, will appear here.
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}

\examples{

cat("See examples in nlmrt-package.Rd\n")


}
\keyword{ nonlinear least squares }

