\name{modgr}
\alias{modgr}
\title{Compute gradient from residuals and Jacobian.
}
\description{
   For a nonlinear model originally expressed as an expression of the form
         lhs ~ formula_for_rhs
   assume we have a resfn and jacfn that compute the residuals and the 
   Jacobian at a set of parameters. This routine computes the gradient, 
   that is, t(Jacobian) . residuals. 
}
\usage{
   modgr(prm, resfn, jacfn, ...)
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
  \item{jacfn}{
     A function to compute the Jacobian of the residuals at a paramter vector.
}
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{modgr} calls resfn to compute residuals and jacfn to compute the
    Jacobian at the parameters \code{prm} using external data in the dot arguments.
    It then computes the gradient using t(Jacobian) . residuals.

    Note that it appears awkward to use this function in calls to optimization
    routines. The author would like to learn why. 
}
\value{
    The numeric vector with the gradient of the sum of squares at the paramters.
}
\references{

     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

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
  y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 
      50.156, 62.948, 75.995, 91.972)  # for testing
  tt <- seq_along(y)  # for testing
  f <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt))
  p <- c(b1 = 1, b2 = 1, b3 = 1)
  myres <- model2resfun(f, p)
  myjac <- model2jacfun(f, p)
  mygr  <- model2grfun(f, p)
  gr <- mygr(p, tt = tt, y = y)
  grm <- modgr(p, myres, myjac, tt = tt, y = y)
  cat("max(abs(grm - gr)) =",max(abs(grm-gr)),"\n")
}
\keyword{ nonlinear least squares }

