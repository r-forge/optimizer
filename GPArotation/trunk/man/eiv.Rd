\name{eiv}
\alias{eiv}

\title{Errors-in-Variables Rotation}
\usage{
    eiv(L, identity=seq(ncol(L)), normalize=FALSE)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{identity}{indicates rows which should be identity matrix.}
    \item{normalize}{not sure if this makes sense}
}

\value{A list as needed by \code{factanal} with elements
   \item{loadings}{Lh from \code{GPForth} or \code{GPFoblq}.} 
   \item{Th}{The rotation} 
   \item{method}{A string indicating the rotation objective function ("eiv").}
   \item{orthogonal}{For consistency with other rotation results. Always FALSE.}
   \item{convergence}{For consistency with other rotation results. Always TRUE.}
   \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.}
   }
\description{
   Rotate to errors-in-variables representation.
}
\details{
 This function rotates to an errors-in-variables representation. The
 optimization is not iterative and does not use the GPA algorithm.  
 The function can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.

 The loadings matrix is rotated so the \eqn{k}{k} rows indicated by \code{identity} 
 form an identity matrix, and the remaining \eqn{M-k}{M-k} rows are free parameters.
 \eqn{\Phi}{Phi} is also free. The default makes the first \eqn{k}{k} rows 
 the identity. If inverting the matrix
 of the rows indicated by \code{identity} fails, the rotation will fail and the
 user needs to supply a different choice of rows.

 Not all authors consider this representation to be a rotation. 
 Viewed as a rotation method, it is oblique, with an
 explicit solution: given an initial loadings matrix \eqn{L}{L} partitioned as
 \eqn{L = (L_1^T, L_2^T)^T}{L = rbind(L1, L2)}, then (for the default 
 \code{identity}) the new loadings matrix is 
 \eqn{(I,  (L_2 L_1^{-1})^T)^T}{rbind(I,  L2 \%*\% solve(L1))} 
 and \eqn{\Phi = L_1 L_1^T}{Phi = L1 \%*\% t(L1)}, where \eqn{I}{I} is the \eqn{k}{k} 
 by \eqn{k}{k} identity matrix. It is
 assumed that \eqn{\Phi = I}{Phi = I} for the initial loadings matrix.

 One use of this parameterization is
 for obtaining good starting values (so it looks a little strange
 to rotate towards this solution afterwards). It has a few other purposes: 
 (1) It can be useful for comparison with
 published results in this parameterization; 
 (2) The
 S.E.s are more straightfoward to compute, because it is the solution 
 to an unconstrained
 optimization (though not necessarily computed as such); 
 (3) One
 may have an idea about which reference variables load on only one
 factor, but not impose restrictive constraints on the other loadings,
 so, in a nonrestrictive
 way, it has similarities to CFA; 
 (4) For some purposes, only the subspace spanned by the factors
 is important, not the specific parameterization within this subspace;
 (5) The back-predicted indicators (explained portion of the indicators) 
 do not depend
 on the rotation method. Combined with the greater ease to obtain
 correct standard errors of this method, this allows easier and more
 accurate prediction-standard errors.  
}

\examples{
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, rotation="none")

  fa.eiv <- eiv(fa.unrotated$loadings)
 
  fa.eiv2 <- factanal(factors = 2, covmat=NetherlandsTV, 
                rotation="eiv", control=list(rotate=list(normalize=TRUE)))
  
  cbind(loadings(fa.unrotated), loadings(fa.eiv), loadings(fa.eiv2))

  fa.eiv3 <- eiv(fa.unrotated$loadings, identity=6:7)
  cbind(loadings(fa.unrotated), loadings(fa.eiv), loadings(fa.eiv3))

  }
\seealso{
   \code{\link{rotations}}, 
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}}
}
\references{
   Yasuo Amemiya and his associates .?

   Hagglund and others??

   Tom Wansbeek and Erik Meijer (2000) \emph{Measurement Error and 
   Latent Variables in Econometrics}, Amsterdam: North-Holland.
  }

\author{Erik Meijer and Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


