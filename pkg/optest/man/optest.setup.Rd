\name{optest.setup}
\alias{optest.setup}
\encoding{UTF-8}
\title{Setup general-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  Setup general-purpose optimization wrapper.
}
\usage{
optest.setup(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=NULL, itnmax=NULL, hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{method}{A list of the methods to be used. 
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}.
       The default of NULL causes an appropriate set of methods to be supplied
       depending on the presence or absence of bounds on the parameters. The default
       unconstrained set is Rvmminu, Rcgminu, lbfgsb3, newuoa and nmkb.
       The default bounds constrained set is Rvmminb, Rcgminb, lbfgsb3, bobyqa and nmkb.}
 \item{itnmax}{If provided as a vector of the same length as the list of methods \code{method}, 
	gives the maximum number of iterations or function values for the corresponding 
	method. If a single number is provided, this will be used for all methods. Note that
	there may be control list elements with similar functions, but this should be the
	preferred approach when using \code{optimx}.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim().}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.


}
\value{
   A structure with the optimization setup.
}
\note{
  Most methods in \code{optimx} will work with one-dimensional \code{par}s, but such
  use is NOT recommended. Use \code{\link{optimize}} or other one-dimensional methods instead.

  There are a series of demos available. Once the package is loaded (via \code{require(optimx)} or
  \code{library(optimx)}, you may see available demos via 

  demo(package="optimx")

  The demo 'brown_test' may be run with the command
  demo(brown_test, package="optimx")


Knowledgeable users may take advantage of these functions if they are carrying
out production calculations where the setup and checks could be run once.

}
\source{

See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

}
\references{

 See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

}

\keyword{nonlinear}
\keyword{optimize}
