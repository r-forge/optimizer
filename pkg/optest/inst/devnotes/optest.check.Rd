\name{optest.check}
\alias{optest.check}
\encoding{UTF-8}
\title{check function and parameters for General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  Check for Test General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optimx} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimx() does not call the SANN method.
}
\usage{
optest.check(par, ufn, ugr=NULL, uhess=NULL, lower=-Inf, upper=Inf, 
            ctrl=NULL, ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{ufn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{ugr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.}
 \item{uhess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{ctrl}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

}
\value{



   If there are \code{npar} parameters, then the result is a dataframe having one row
   for each method for which results are reported, using the method as the row name,
   with columns

   \code{par_1, .., par_npar, value, fevals, gevals, niter, convcode, kkt1, kkt2, xtimes}

  where
  \describe{
  \item{par_1}{ .. }
  \item{par_npar}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{fevals}{The number of calls to \code{fn}.}
  \item{gevals}{The number of calls to \code{gr}. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.}
  \item{niter}{For those methods where it is reported, the number of ``iterations''. See
    the documentation or code for particular methods for the meaning of such counts.}
  \item{convcode}{An integer code. \code{0} indicates successful
    convergence. Various methods may or may not return sufficient information
	to allow all the codes to be specified. An incomplete list of codes includes
    \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
      \item{\code{10}}{indicates degeneracy of the Nelder--Mead simplex.}
      \item{\code{51}}{indicates a warning from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{52}}{indicates an error from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
    }
  }
  \item{kkt1}{A logical value returned TRUE if the solution reported has a ``small'' gradient.}
  \item{kkt2}{A logical value returned TRUE if the solution reported appears to have a 
  positive-definite Hessian.}
  \item{xtimes}{The reported execution time of the calculations for the particular method.}
  }


The attribute "details" to the returned answer object contains information,
if computed, on the gradient (\code{ngatend}) and Hessian matrix (\code{nhatend}) 
at the supposed optimum, along with the eigenvalues of the Hessian (\code{hev}), 
as well as the \code{message}, if any, returned by the computation for each \code{method},
which is included for each row of the \code{details}. 
If the returned object from optimx() is \code{ans}, this is accessed 
via the construct
    \code{attr(ans, "details")}

This object is a  matrix based on a list so that if ans is the output of optimx
then attr(ans, "details")[1, ] gives the first row and 
attr(ans,"details")["Nelder-Mead", ] gives the Nelder-Mead row. There is 
one row for each method that has been successful 
or that has been forcibly saved by save.failures=TRUE. 

There are also attributes
   \describe{
   \item{maximize}{to indicate we have been maximizing the objective}
   \item{npar}{to provide the number of parameters, thereby facilitating easy
        extraction of the parameters from the results data frame}
   \item{follow.on}{to indicate that the results have been computed sequentially,
        using the order provided by the user, with the best parameters from one
        method used to start the next. There is an example (\code{ans9}) in 
        the script \code{ox.R} in the demo directory of the package.}
   }
}
\note{
  Most methods in \code{optimx} will work with one-dimensional \code{par}s, but such
  use is NOT recommended. Use \code{\link{optimize}} or other one-dimensional methods instead.

  There are a series of demos available. Once the package is loaded (via \code{require(optimx)} or
  \code{library(optimx)}, you may see available demos via 

  demo(package="optimx")

  The demo 'brown_test' may be run with the command
  demo(brown_test, package="optimx")

  The package source contains several functions that are not exported in the
  NAMESPACE. These are 
  \describe{
  \item{\code{optimx.setup()}}{ which establishes the controls for a given run;}
  \item{\code{optimx.check()}}{ which performs bounds and gradient checks on
      the supplied parameters and functions;}
  \item{\code{optimx.run()}}{which actually performs the optimization and post-solution
      computations;}
  \item{\code{scalecheck()}}{ which actually carries out a check on the relative scaling
      of the input parameters.}
  }

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

 Nash JC (2014). On Best Practice Optimization Methods in R., 
        \emph{Journal of Statistical Software}, 60(2), 1-14.,
        URL http://www.jstatsoft.org/v60/i02/.

}

\keyword{nonlinear}
\keyword{optimize}
