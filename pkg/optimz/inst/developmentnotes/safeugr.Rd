\name{ugr}
\alias{ugr}
\encoding{UTF-8}
\title{Wrapper for user gradient function for optimization tools}
\concept{minimization}
\concept{maximization}
\description{
        Provides a wrapper around user gradient function for nonlinear optimization
	to try to control for inadmissible arguments to user objective, gradient
	or hessian functions, as well as provide for maximization.
}
\usage{
ugr(par, fnuser)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
}
\details{
   The usual dot arguments (...) are subsumed in fnuser$dots to save
   complexity in the function call. Note that we need to unlist()
   these is the call to the actual user function.
}
\value{
  \code{ugr} returns a vector numeric value, but all elements are set to the R 
   constant .Machine$double.xmax if the inputs to the function are inadmissible and the
   computation of \code{gr} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}
\examples{
cat("All examples in tests/ subdirectory\n")

}

\keyword{nonlinear}
\keyword{optimize}
