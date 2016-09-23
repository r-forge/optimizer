% File src/library/stats/man/optim.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optimr}
\alias{optimr}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optim} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimrx::optimr() does not call the SANN method.
}
\usage{
optimr(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE,
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

   If 'gr' is \code{NULL}, whatever default action is specified for the 
   chosen method for the case of a null gradient code
   is used. For many methods, this is a finite-difference approximation, 
   but some methods require user input for the gradient. 
   
   If 'gr' is a character string, then that string is taken as the name of 
   a gradient approximation function, for example, "grfwd", "grback" and
   "grcentral" for standard forward, backward and central approximations.
   (See the code in package \code{optextras} for details.) Method "grnd" 
   uses the \code{grad()} function from package \code{numDeriv}.

 }   
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints. A small set of methods can handle masks, that is,
   fixed parameters, and these can be specified by making the lower and upper bounds
   equal to the starting value. (It is possible that the starting value could be different
   from the lower/upper bounds set,
   but this behaviour has NOT yet been defined and users are cautioned.)}
 \item{method}{A character string giving the name of the optimization method to be
    applied.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. Note that this will NOT necessarily
       use the same approximation as may be provided by the method called. Instead, 
       the function \code{hessian()} from package \code{numDeriv} is used if no gradient
       \code{gr} is supplied, else the function \code{jacobian()} from \code{numDeriv}
       is applied to the gradient function \code{gr}.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{Further arguments to be passed to \code{fn} 
    and \code{gr} if needed for computation of these quantities; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} should be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this. DO NOT
  use both mechanisms simultaneously. 

  Possible method codes are 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 
  'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 
  'lbfgsb3', 'Rtnmin', 'hjn', 'lbfgs' or 'Rvmmin'. These are in base R or in CRAN repositories. 

  If no method is specified, "Nelder-Mead" will be attempted.
  
  Function \code{fn} must return a finite value at the initial set
  of parameters. Some methods can handle \code{NA} or \code{Inf} 
  if the function cannot be evaluated at the supplied value. However, some methods, of
  which \code{"L-BFGS-B"} is known to be a case, require that the values
  returned should always be finite.

  While \code{optim} can be used recursively, and for a single parameter
  as well as many, this may not be true for other methods in \code{optimrx}. 
  \code{optim} also accepts a zero-length \code{par}, and just evaluates the function 
  with that argument. 

  Generally, you are on your own if you choose to apply constructs mentioned in the
  above two paragraphs.

  For details of methods, please consult the documentation of the individual methods.
  (The NAMESPACE file lists the packages from which functions are imported.)
  However, method \code{"hjn"} is a conservative implementation of a Hooke and 
  Jeeves (1961) and is part of this package.

  The \code{control} argument is a list that can supply any of the
  following components (and may supply others eventually):

  \describe{
    \item{\code{trace}}{Non-negative integer. If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing. trace = 0 gives no output 
      (To understand exactly what these do see the source code: higher 
      levels give more detail.)}
    \item{\code{save.failures}}{ = TRUE (default) if we wish to keep "answers" from runs 
      where the method does not return convcode==0. FALSE otherwise.}
    \item{\code{maximize}}{ = TRUE if we want to maximize rather than minimize 
      a function. (Default FALSE). }
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods. Default FALSE.}
    \item{\code{dowarn}}{= TRUE if we want warnings generated by optimx. Default is 
      TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.
      Default is (0.5)*.Machine$double.xmax }
  }

  NOte that some \code{control} elements apply only to some of methods. 
  See individual packages for details. 

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Apparently no other attributes of \code{par}
  are copied over. (We have not verified this as at 2009-07-29.)

}
\value{

    A list with components:
  \describe{
   \item{par}{The best set of parameters found.}
   \item{value}{The value of ‘fn’ corresponding to ‘par’.}
   \item{counts}{ A two-element integer vector giving the number of calls to
          ‘fn’ and ‘gr’ respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to ‘fn’
          to compute a finite-difference approximation to the gradient.}
   \item{convergence}{ An integer code. ‘0’ indicates successful completion. }
   \item{ message}{ A character string giving any additional information returned
          by the optimizer, or ‘NULL’.}
   \item{hessian}{If requested, an approximation to the hessian of ‘fn’
        at the final parameters.}
  }
  }
\source{
See the manual pages for \code{optim()} and the packages listed in NAMESPACE.
}
\references{

See the manual pages for \code{optim()} and the packages listed in NAMESPACE.

 Hooke R. and Jeeves, TA (1961). Direct search solution of numerical and statistical problems. 
   Journal of the Association for Computing Machinery (ACM). 8 (2): 212–229.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

}
\keyword{nonlinear}
\keyword{optimize}
