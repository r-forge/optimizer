\name{uhess}
\alias{uhess}
\encoding{UTF-8}
\title{Wrapper for user Hessian function for optimization tools}
\concept{minimization}
\concept{maximization}
\description{
Provides a wrapper around user analytic Hessian function for nonlinear optimization
	to try to control for inadmissible arguments to user that function,
	as well as provide for maximization.
}
\usage{
uhess(par, fnuser)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
}
\details{ Note that a FIXED name environment opx12env is used to carry the 
   user's fn, gr and hess and the counters.
}
\value{
  \code{uhess} returns a matrix numeric value, but if the inputs to the function
    are inadmissible, all elements are set to the R 
   constant .Machine$double.xmax and the
   computation of \code{fn} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}\examples{
cat("uhess: try an inadmissible set of parameters to a user function\n")

}
\keyword{nonlinear}
\keyword{optimize}
