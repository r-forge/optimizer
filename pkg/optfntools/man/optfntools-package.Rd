\name{optfntools-package}
\alias{optfntools-package}
\alias{optfntools}
\docType{package}
\title{
A set of functions that wrap user objective, gradient and hessian functions, that
allow the user functions to be checked for consistency and validity, that parameters
satisfy bounds, or which check a set of parameters for optimality conditions, including
generation of an appropriate numerical approximation to gradient or hessian.
}
\description{
The three functions ufn, ugr and uhess wrap corresponding user functions fn, gr, and 
hess so that these functions can be executed safely (via try()) and also so parameter or
function scaling can be applied. The wrapper
functions also allow for maximization of functions (via minimization of the negative of
the function) using the logical parameter \code{maximize}.

There are three test functions, fnchk, grchk, and hesschk, to allow the user 
function to be tested for validity and correctness. However, no set of tests is 
exhaustive, and extensions and improvements are welcome. The package 
\code{numDeriv} is used for generation of numerical approximations to 
derivatives.
}
\details{
\tabular{ll}{?? fix ??\cr
Package: \tab optfntools\cr
Version: \tab 2011-9.25\cr
Date: \tab 2011-09-25\cr
License: \tab GPL-2\cr
Lazyload: \tab Yes\cr
Repository: \tab R-Forge\cr
Repository/R-Forge/Project: \tab optimizer\cr
Repository/R-Forge/Revision: \cr
Date/Publication: \cr
}

Index:
\preformatted{
ufn          Wrapper for user objective function
ugr          Wrapper for user gradient function
uhess        Wrapper for user Hessian function
fnchk        Test validity of user function
grchk        Check that gradient function evaluation 
             matches numerical gradient
hesschk      Check that Hessian function evaluation 
             matches numerical approximation
gHgen        Compute gradient and Hessian as a given 
             set of parameters
gHgenb       Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
ugHgenb      Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
             but using the fnuser list of fn, gr, and 
             hess.
}
}
\author{
John C Nash <nashjc@uottawa.ca>

Maintainer: John C Nash <nashjc@uottawa.ca>
}

\keyword{ optimization }
\seealso{
optimx
}

