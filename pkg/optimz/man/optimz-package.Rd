\name{optimz-package}
\alias{optimz-package}
\docType{package}
\title{
An enlarged replacement and extension of the optim() function
}
\description{

\code{optimz}

provides a replacement and extension of the \code{link{optim()}} function to 
unify and streamline optimization capabilities in R for smooth, possibly box
constrained functions of several or many parameters. The \code{optimz} package
also subsumes the packages \code{optextras}, \code{Rcgmin}, \code{Rvmmin}
 and \code{Rtnmin}.

Examples can be found in the \code{tests/} directory, and in the 
demo directory, \code{demo(package = "optimz")}. ??

The three functions \code{ufn}, \code{ugr} and \code{uhess}
wrap corresponding user functions \code{fn}, \code{gr}, and 
\code{hess} so that these functions can be executed safely (via 
\code{try()}) and also so parameter or
function scaling could be applied. However at 2015-4-27, automated scaling has been found to
be very time consuming and has been left out. The wrapper
functions also allow for maximization of functions (via minimization of the negative of
the function) using the logical parameter \code{maximize}.

There are three test functions, fnchk, grchk, and hesschk, to allow the user 
function to be tested for validity and correctness. However, no set of tests is 
exhaustive, and extensions and improvements are welcome. The package 
\code{numDeriv} is used for generation of numerical approximations to 
derivatives.

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
\details{

Index:
\preformatted{
axsearch     Perform an axial search optimality check
bmchk        Check bounds and masks for parameter constraints
bmstep       Compute the maximum step along a search direction.
fnchk        Test validity of user function
gHgen        Compute gradient and Hessian as a given 
             set of parameters
gHgenb       Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
grback       Backward numerical gradient approximation
grcentral    Central numerical gradient approximation
grchk        Check that gradient function evaluation 
             matches numerical gradient
grfwd        Forward numerical gradient approximation
grnd         Gradient approximation using \code{numDeriv}
hesschk      Check that Hessian function evaluation 
             matches numerical approximation
kktc         Check the Karush-Kuhn-Tucker optimality conditions
scalecheck   Check scale of initial parameters and bounds
ufn          Wrapper for user objective function
ugr          Wrapper for user gradient function
ugHgenb      Compute gradient and Hessian as a given 
             set of parameters appying bounds and masks
             but using the opx12env list of fn, gr, and 
             hess. Note FIXED name opx12env
uhess        Wrapper for user Hessian function
}
}
\keyword{ package }
\keyword{ optimization }
\seealso{optim}

