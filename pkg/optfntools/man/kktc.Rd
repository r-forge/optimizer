% File src/library/stats/man/kktc.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{kktc}
\alias{kktc}
\encoding{UTF-8}
\title{Check Kuhn Karush Tucker conditions for a supposed function minimum}
\concept{minimization}
\concept{maximization}
\description{
  Provide a check on Kuhn-Karush-Tucker conditions based on quantities
  already computed. Some of these used only for reporting.
}
\usage{
     kktc(par, fval, ngr, nHes, nbm, maxfn=FALSE, control=list() )
}
\arguments{
 \item{par}{A vector of values for the parameters which are supposedly optimal.}
 \item{fval}{Function value at parameters \code{par}.}
 \item{ngr}{Gradient (possibly estimated) of the function at \code{par}.}
 \item{nHes}{Matrix that is supposedly the Hessian at \code{par}.}
 \item{nbm}{Number of active masks or bounds on \code{par}.}
 \item{maxfn}{Logical TRUE if function is being maximized. Default FALSE.}
 \item{control}{A list of controls for the function, made up of the following items}
    \itemize{
       \item{kkttol}{A tolerance for the KKT gradient check.}
       \item{kkt2tol}{A tolerance used for the KKT Hessian check.}
       \item{ktrace}{A logical flag that, if TRUE, displays intermediate information
          on the progress of the function. Default is FALSE.}
    }
}
\details{
  At the moment, kktc does NOT handle box constraints or masks in its
  tests, except to ignore the dimensions for masks (not active bounds).
  The present code was developed by lifting it from \code{optimx}
  so that KKT checks could be added to other optimization routines.
  Note that we IGNORE the sign of the gradient in the case that the
  objective is being maximized. We do account for the sign of the
  Hessian for curvature however.
}

\value{
  The output is a list consisting of 
  \item{gmax}{The absolute value of the largest gradient component in magnitude.}
  \item{evratio}{The ratio of the smallest to largest Hessian eigenvalue. Note that this
     may be negative.}
  \item{kkt1}{A logical value that is TRUE if we consider the first (i.e., gradient) 
     KKT condition to be satisfied. WARNING: The decision is dependent on tolerances and
     scaling that may be inappropriate for some problems.}
  \item{kkt2}{A logical value that is TRUE if we consider the second (i.e., positive
      definite Hessian) KKT condition to be satisfied. WARNING: The decision is dependent 
      on tolerances and scaling that may be inappropriate for some problems.}
}
\seealso{

    \code{\link{optim}}

}
\examples{
# genrose function code
}
\keyword{nonlinear}
\keyword{optimize}
