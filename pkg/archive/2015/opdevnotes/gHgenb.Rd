% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{gHgenb}
\alias{gHgenb}
\encoding{UTF-8}
\title{Generate gradient and Hessian for a function at given parameters.}
\description{
  \code{gHgenb} is used to generate the gradient and Hessian of an objective
  function used for optimization. If a user-provided gradient function 
  \code{gr} is available it is used to compute the gradient, otherwise 
  package \code{numDeriv} is used. If a user-provided Hessian function
  \code{hess} is available, it is used to compute a Hessian. Otherwise, if
  \code{gr} is available, we use the function \code{jacobian()} from
  package \code{numDeriv} to compute the Hessian. In both these cases we
  check for symmetry of the Hessian. Computational Hessians are commonly
  NOT symmetric. If only the objective function \code{fn} is provided, then
  the Hessian is approximated with the function \code{hessian} from 
  package \code{numDeriv} which guarantees a symmetric matrix. 
}
\usage{
  gHgenb(par, fn, gr=NULL, hess=NULL, bdmsk=NULL, lower=NULL, upper=NULL,
      control=list(ktrace=0), ...) 
}
\arguments{
 \item{par}{Set of parameters, assumed to be at a minimum of the function \code{fn}.}
 \item{fn}{Name of the objective function.}
 \item{gr}{(Optional) function to compute the gradient of the objective function. If present,
     we use the Jacobian of the gradient as the Hessian and avoid one layer of numerical
     approximation to the Hessian.} 
 \item{hess}{(Optional) function to compute the Hessian of the objective function. This
     is rarely available, but is included for completeness.}
 \item{bdmsk}{An integer vector of the same length as \code{par}. When an element
     of this vector is 0, the corresponding parameter value is fixed (masked) 
     during an optimization. Non-zero values indicate a parameter is free (1),
     at a lower bound (-3) or at an upper bound (-1), but this routine only
     uses 0 values.}
 \item{lower}{Lower bounds for parameters in \code{par}.}
 \item{upper}{Upper bounds for parameters in \code{par}.}
 \item{control}{A list of controls to the function. Currently 
                  asymptol (default of 1.0e-7 which tests for asymmetry of Hessian approximation
                  (see code for details of the test); 
                  ktrace, a logical flag which, if TRUE, monitors the progress 
                  of gHgenb (default FALSE), and 
                  stoponerror, defaulting to FALSE to NOT stop when there is
                  an error or asymmetry of Hessian. Set TRUE to stop.}
 \item{\dots}{Extra data needed to compute the function, gradient and Hessian.}
}

\details{
  None
}
\value{
  \code{ansout} a list of four items, 
  \itemize{
      \item{\code{gn}}{  The approximation to the gradient vector.}
      \item{\code{Hn}}{  The approximation to the Hessian matrix.}
      \item{\code{gradOK}}{  TRUE if the gradient has been computed acceptably. FALSE otherwise.}
      \item{\code{hessOK}}{  TRUE if the gradient has been computed acceptably and passes the
             symmetry test. FALSE otherwise.}
      \item{\code{nbm}}{  The number of active bounds and masks.}
  }
}
\examples{
cat("all examples in tests/ subdirectory\n")

}
\keyword{nonlinear}
\keyword{optimize}
