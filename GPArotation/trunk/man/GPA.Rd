\name{GPA}
\alias{GPForth}
\alias{GPFoblq}

\title{Rotation optimization}
\usage{
    GPForth(A, Tmat=diag(ncol(A)), method="varimax", normalize=FALSE, eps=1e-5, maxit=1000, ...)
    GPFoblq(A, Tmat=diag(ncol(A)), method="quartimin", normalize=FALSE, eps=1e-5, maxit=1000, ...)
    }
\arguments{
    \item{A}{initial factor loadings matrix for which the rotation 
             criterian is to be optimized.}
    \item{Tmat}{initial rotation matrix.}
    \item{method}{rotation objective criterian.}
    \item{normalize}{see details.}
    \item{eps}{convergence is assumed when the norm of the gradient is smaller
         than eps.}
    \item{maxit}{maximum number of iterations allowed in the main loop.}
    \item{...}{additional arguments passed to the rotation objective}
}

\description{
   Gradient projection rotation optimization routine used by various rotation
   objective.
   }
\value{A GPArotation object which is a list with elements 
  \item{Lh}{The rotated loadings, one column for each factor.}
  \item{Th}{The rotation matrix,  Lh \%*\% t(Th) = A.}
  \item{Table}{A matrix recording the iterations of the rotation optimization.}
  \item{method}{A string indicating the rotation objective function.}
  \item{orthogonal}{A logical indicating if the rotation is orthogonal.}
  \item{convergence}{A logical indicating if convergence was obtained.}
  \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.
  This will be the identity matrix for orthogonal
  rotations so is omitted from the result of GPForth.}
  }
\details{
 Gradient projection rotation optimization routines developed by
 Coen A. Bernaards and Robert I. Jennrich. These functions can be used directly
 to rotate a loadings matrix, or indirectly through a rotation
 objective passed to a factor estimation routine such as \code{\link{factanal}}.
 For examples of the indirect use see the documention for 
 rotations (such as \code{\link{oblimin}}).

 GPForth is the main GP algorithm for orthogonal rotation.
 GPFoblq is the main GP algorithm for oblique rotation.
 Both algorithms require a loadings matrix \code{A} which fixes the 
 equivalence class over which the optimization is done. It must be 
 the solution to the orthogonal factor analysis problem. 
 A rotation is defined as code{A \%*\% t(solve(Tmat))}. For the 
 orthogonal case \code{Tmat} is orthonormal so this simplifies 
 to code{A \%*\% Tmat}. The starting point for 
 iterative optimization is given by the \code{Tmat} rotation of \code{A}. 
 By default the initial rotation is the 
 identity matrix. For some rotation criteria local optima may exist
 and it is recommended to check for this by starting with many different
 initial rotations. The function \code{\link{Random.Start}} will help to do
 this.
 
 The argument \code{method} can be used to specify a string indicating
 the rotation objective. \code{GPFoblq} defaults to \code{"quartimin"} 
 and \code{GPForth} defaults to \code{"varimax"}. Available rotation objectives
 are "oblimin", "quartimin", "target", "pst", "oblimax", "entropy", 
 "quartimax", "varimax", "simplimax", "bentler", "tandemI", "tandemII", 
 "geomin", "cf", "infomax" and "mccammon". The string is prefixed with 
 "vgQ." to give the actual function call. The \code{vgQ.*} function call
 would typically not be used directly, so these methods are not exported
 from the package namespace. You can print these functions to see the code for
 calculating a criterion, but since they are not exported the package 
 name needs to be specified. For example, use
 \code{GPArotation:::vgQ.oblimin} to view the function \code{vgQ.oblimin}.

 Some rotation criteria (including \code{simplimax}, \code{pst}, 
 \code{procrustes})
 require one or more additional arguments. For example, \code{simplimax}
 needs the number of 'close to zero loadings' which is given  as
 the extra argument \code{k}. Check the rotation methods for details. 
 (If a new rotation method is implemented and needs additional
 arguments then this is the  way to pass them.)

 The argument \code{normalize} gives an indication of if and how any
 normalization should be done before rotation, and then undone after rotation.
 If \code{normalize} is \code{FALSE} (the default) no normalization is done.
 If \code{normalize} is \code{TRUE} then Kaiser normalization is done. (So
 squared row entries of normalized \code{A} sum to 1.0.)
 If \code{normalize} is a vector of length equal to the number of indicators (=
 number of rows of \code{A}) then the colums are divided by \code{normalize} 
 before rotation and multiplied by \code{normalize} after rotation. 
 If \code{normalize} is a function then it should take \code{A} as an argument
 and return a vector which is used like the vector above. 
}
\seealso{
   \code{\link{Random.Start}}
   \code{\link[stats]{factanal}}
   \code{\link{oblimin}}
   \code{\link{quartimin}}
   \code{\link{targetT}}
   \code{\link{targetQ}}
   \code{\link{pstT}}
   \code{\link{pstQ}}
   \code{\link{oblimax}}
   \code{\link{entropy}}
   \code{\link{quartimax}}
   \code{\link{Varimax}}
   \code{\link[stats]{varimax}}
   \code{\link{simplimax}}
   \code{\link{bentlerT}}
   \code{\link{bentlerQ}}
   \code{\link{tandemI}}
   \code{\link{tandemII}}
   \code{\link{geominT}}
   \code{\link{geominQ}}
   \code{\link{cfT}}
   \code{\link{cfQ}}
   \code{\link{infomaxT}}
   \code{\link{infomaxQ}}
   \code{\link{mccammon}}
   \code{\link[stats]{promax}}
}
\examples{
  data("Harman", package="GPArotation")
  qHarman  <- GPForth(Harman8, Tmat=diag(2), method="quartimax")

  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, 
              normalize=TRUE, rotation="none")

  GPForth(loadings(fa.unrotated), method="varimax", normalize=TRUE)$Lh

  TV <- GPFoblq(loadings(fa.unrotated), method="oblimin", normalize=TRUE)

  print(TV)
  print(TV, Table=TRUE)
  summary(TV)
  }
\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert}
\source{code is modified from original source splusfunctions.net available at 
  \url{http://www.stat.ucla.edu/research/gpa}.
}
\reference{
   source{\url{http://www.stat.ucla.edu/research} or
          \url{http://www.stat.ucla.edu/research/gpa}
   }
   The software reference is

   Bernaards, C.A. and Jennrich, R.I. (in press) Gradient Projection 
     Algorithms and Software for Arbitrary Rotation Criteria in Factor
     Analysis. Educational and Psychological Measurement.

   Theory of gradient projection algorithms may be found in:

   Jennrich, R.I. (2001). A simple general procedure for orthogonal rotation.
       Psychometrika, 66, 289-306.

   Jennrich, R.I. (2002). A simple general method for oblique rotation. 
       Psychometrika, 67, 7-19.
   }
\concept{rotation}
\keyword{multivariate}
