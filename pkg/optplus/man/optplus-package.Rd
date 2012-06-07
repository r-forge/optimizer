\name{optplus-package}
\alias{optplus-package}
\docType{package}
\title{
A replacement and extension of the optim() function, plus various optimization
tools
}
\description{
Provides a replacement and extension of the optim() function to unify
and streamline optimization capabilities in R for smooth, possibly box
constrained functions of several or many parameters

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
\tabular{ll}{
Package: \tab optplus\cr
Version: \tab 2012-6.07\cr
Date: \tab 2012-06-07\cr
License: \tab GPL-2\cr
Lazyload: \tab Yes\cr
Depends: \tab numDeriv \cr
Suggests: \tab BB, ucminf, Rcgmin, Rvmmin, minqa, setRNG, dfoptim\cr
Repository: \tab R-Forge\cr
Repository/R-Forge/Project: \tab optimizer\cr
}

Index:
\preformatted{
axsearch     Perform an axial search optimality check
bmchk        ???
bmstep       ???
fnchk        Test validity of user function
get.best     ???
get.result   Result extraction tool for answers from optimx
gHgen        Compute gradient and Hessian as a given 
             set of parameters
gHgenb       Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
grback       ???
grcentral    ???
grchk        Check that gradient function evaluation 
             matches numerical gradient
grfwd        ???
grnd         ???
hesschk      Check that Hessian function evaluation 
             matches numerical approximation
kktc         Check the Karush-Kuhn-Tucker optimality conditions
optansout    Output tool for answers from optimx
optimx       General-purpose optimization wrapper
scalecheck   ??
ufn          Wrapper for user objective function
ugr          Wrapper for user gradient function
ugHgenb      Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
             but using the opx12env list of fn, gr, and 
             hess. Note FIXED name opx12env
uhess        Wrapper for user Hessian function
}
}
\author{
John C Nash <nashjc@uottawa.ca> and Ravi Varadhan <RVaradhan@jhmi.edu>

Maintainer: John C Nash <nashjc@uottawa.ca>
}
\references{
  Nash, John C. and Varadhan, Ravi (2011) Unifying Optimization Algorithms 
          to Aid Software System Users: optimx for R, Journal of Statistical
	  Software, publication pending.
}
\keyword{ package }
\keyword{ optimization }
\seealso{
optim
}

