\name{cholesky}
\alias{cholesky}

\title{Cholesky Rotation}
\usage{
    cholesky(L, Phi=NULL, reference=seq(NCOL(L)), ...)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{Phi}{factor covariance matrix corresponding to L}
    \item{reference}{indicates rows of loading matrix that should be 
       used as reference values (to which Cholesky decompositon is applied).}
    \item{...}{additional arguments discarded.}
}

\value{A list (which includes elements used by \code{factanal}) with:
   \item{loadings}{The new loadings matrix.} 
   \item{Th}{The rotation.} 
   \item{method}{A string indicating the rotation objective function ("cholesky").}
   \item{orthogonal}{For consistency with other rotation results. 
     This is set to TRUE if \code{Phi} is NULL (default) and FALSE otherwise.}
   \item{convergence}{For consistency with other rotation results. Always TRUE.}
   \item{Phi}{The covariance matrix of the rotated factors. 
      NULL for orthogonal solutions.}
   }
\description{
   Rotate using cholesky parameterization.
}
\details{
 This function rotates to a Cholesky representation. The
 optimization is not iterative and does not use the GPA algorithm.  
 The function can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.

 The loadings matrix is rotated so the \eqn{k}{k} rows (indicator variables)
 indicated by \code{reference} have a Cholesky factorization, and the remaining 
 \eqn{M-k}{M-k} rows are rotated as necessary to give a loadings matrix 
 equivalent to \code{L}.

 If the argument \eqn{\Phi}{Phi} is absent or null, it is assumed that it
 is the identity (orthogonal solution). In this case the flag \code{orthogonal}
 in the result is set to \code{TRUE} and otherwise it is set to \code{FALSE}
 however, IS THE NEXT CORRECT? strictly speaking, the result is orthogonal if
 the original is, and not if the original is not.

 The default uses the first \eqn{k}{k} rows 
 the as the reference. If the Cholesky decomposition for the matrix
 of the rows indicated by \code{reference} fails, the rotation will fail and the
 user needs to supply a different choice of rows.

 One use of this parameterization is
 for obtaining good starting values (so it may appear strange
 to rotate towards this solution afterwards). It has a few other purposes: 
 (1) It can be useful for comparison with
 published results in this parameterization; 
 (2) The S.E.s are more straightfoward to compute, because it is the solution 
 to an unconstrained
 optimization (though not necessarily computed as such); 
 (3) The models with k and (k+1) factors are nested, so it
 is more straightforward to test the k-factor model versus the
 (k+1)-factor model. In particular, in addition to the LR test
 (which does not depend on the rotation), now the Wald test
 and LM test can be used as well. For these, the test of a
 k-factor model versus a (k+1)-factor model is a joint test 
 whether all the free parameters (loadings) in the (k+1)st 
 column of \code{L} are zero.
 (4) For some purposes, only the subspace spanned by the factors
 is important, not the specific parameterization within this subspace;
 (5) The back-predicted indicators (explained portion of the indicators) 
 do not depend
 on the rotation method. Combined with the greater ease to obtain
 correct standard errors of this method, this allows easier and more
 accurate prediction-standard errors.  
}

\examples{

#THIS IS EIV EXAMPLE  -- IS THERE SOMETHING BETTER?
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, rotation="none")

  fa.chol <- cholesky(fa.unrotated$loadings)
 
  fa.chol2 <- factanal(factors = 2, covmat=NetherlandsTV, rotation="cholesky")
  
  cbind(loadings(fa.unrotated), loadings(fa.chol), loadings(fa.chol2))

  fa.chol3 <- cholesky(fa.unrotated$loadings, reference=6:7)
  cbind(loadings(fa.unrotated), loadings(fa.chol), loadings(fa.chol3))

  }
\seealso{
   \code{\link{eiv}}, 
   \code{\link{rotations}}, 
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}}
}
\references{
THESE ARE EIV RTEFERENCES  -- THEY ARE PROBABLY NOT RIGHT?
  G{\"{o}}sta H{\"{a}}gglund. (1982). "Factor Analysis by
    Instrumental Variables Methods." \emph{Psychometrika}, 47, 209--222.

   Sock-Cheng Lewin-Koh and Yasuo Amemiya. (2003). "Heteroscedastic factor
   analysis." \emph{Biometrika}, 90, 85--97.

   Tom Wansbeek and Erik Meijer (2000) \emph{Measurement Error and 
   Latent Variables in Econometrics}, Amsterdam: North-Holland.
  }

\author{Erik Meijer and Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


