\name{echelon}
\alias{echelon}

\title{Echelon Rotation}
\usage{
    echelon(L, reference=seq(NCOL(L)), ...)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{reference}{indicates rows of loading matrix that should be 
       used as reference values (to which Cholesky decompositon is applied).}
    \item{...}{additional arguments discarded.}
}

\value{A list (which includes elements used by \code{factanal}) with:
   \item{loadings}{The new loadings matrix.} 
   \item{Th}{The rotation.} 
   \item{method}{A string indicating the rotation objective function ("echelon").}
   \item{orthogonal}{TRUE}
   \item{convergence}{For consistency with other rotation results. Always TRUE.}
   }
\description{
   Rotate using Cholesky parameterization.
}
\details{
 This function rotates to a Cholesky representation. The
 optimization is not iterative and does not use the GPA algorithm.  
 The function can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.
 An orthogonal solution is assumed (so \eqn{\Phi}{Phi} is identity).

 The loadings matrix is rotated so the \eqn{k}{k} rows (indicator variables)
 indicated by \code{reference} have a Cholesky factorization, and the remaining 
 \eqn{M-k}{M-k} rows are rotated as necessary to give a loadings matrix 
 equivalent to \code{L}.

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

 (3) The models with \eqn{k}{k} and \eqn{(k+1)}{(k+1)} factors are nested, so it
 is more straightforward to test the \eqn{k}{k}-factor model versus the
 \eqn{(k+1)}{(k+1)}-factor model. In particular, in addition to the LR test
 (which does not depend on the rotation), now the Wald test
 and LM test can be used as well. For these, the test of a
 \eqn{k}{k}-factor model versus a \eqn{(k+1)}{(k+1)}-factor model is a 
 joint test 
 whether all the free parameters (loadings) in the \eqn{(k+1)}{(k+1)}st 
 column of \code{L} are zero.

 (4) For some purposes, only the subspace spanned by the factors
 is important, not the specific parameterization within this subspace;

 (5) The back-predicted indicators (explained portion of the indicators) 
 do not depend
 on the rotation method. Combined with the greater ease to obtain
 correct standard errors of this method, this allows easier and more
 accurate prediction-standard errors.  

 (6) This parameterization and its standard errors can be used to
 detect identification problems (McDonald, 1999, pp. 181-182).
}

\examples{
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, rotation="none")

  fa.chol <- echelon(fa.unrotated$loadings)
 
  fa.chol2 <- factanal(factors = 2, covmat=NetherlandsTV, rotation="echelon")
  
  cbind(loadings(fa.unrotated), loadings(fa.chol), loadings(fa.chol2))

  fa.chol3 <- echelon(fa.unrotated$loadings, reference=6:7)
  cbind(loadings(fa.unrotated), loadings(fa.chol), loadings(fa.chol3))
  }
\seealso{
   \code{\link{eiv}}, 
   \code{\link{rotations}}, 
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}}
}
\references{
   Roderick P. McDonald (1999)
     \emph{Test Theory: A Unified Treatment}, Mahwah, NJ: Erlbaum

   Tom Wansbeek and Erik Meijer (2000) \emph{Measurement Error and 
   Latent Variables in Econometrics}, Amsterdam: North-Holland.
  }

\author{Erik Meijer and Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


