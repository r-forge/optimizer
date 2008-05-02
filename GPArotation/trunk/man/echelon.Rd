\name{echelon}
\alias{echelon}

\title{Echelon Rotation}
\usage{
    echelon(L, reference=seq(NCOL(L)), ...)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{reference}{indicates rows of loading matrix that should be 
       used to determine the rotation transformation.}
    \item{...}{additional arguments discarded.}
}

\value{A list (which includes elements used by \code{factanal}) with:
   \item{loadings}{The new loadings matrix.} 
   \item{Th}{The rotation.} 
   \item{method}{A string indicating the rotation objective function ("echelon").}
   \item{orthogonal}{For consistency with other rotation results. Always TRUE.}
   \item{convergence}{For consistency with other rotation results. Always TRUE.}
   }
\description{
   Rotate to an echelon parameterization.
}
\details{
 The loadings matrix is rotated so the \eqn{k}{k} rows of the loading matrix 
 indicated by \code{reference} are the Cholesky factorization given by 
 \code{t(chol(L[reference,] \%*\% t(L[reference,])))}.
 This defines the
 rotation transformation, which is then also applied to other rows to 
 give the new loadings matrix.

 The optimization is not iterative and does not use the GPA algorithm.  
 The function can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.
 An orthogonal solution is assumed (so \eqn{\Phi}{Phi} is identity).

 The default uses the first \eqn{k}{k} rows 
 as the reference. If the submatrix of \code{L} indicated by reference is 
 singular then the rotation will fail and the
 user needs to supply a different choice of rows.

 One use of this parameterization is
 for obtaining good starting values (so it may appear strange
 to rotate towards this solution afterwards). It has a few other purposes: 

 (1) It can be useful for comparison with
 published results in this parameterization.

 (2) The S.E.s are more straightforward to compute, because it is the solution 
 to an unconstrained optimization (though not necessarily computed as such). 

 (3) The models with k and (k+1) factors are nested, so it
 is more straightforward to test the k-factor model versus the
 (k+1)-factor model. In particular, in addition to the LR test
 (which does not depend on the rotation), now the Wald test
 and LM test can be used as well. For these, the test of a
 k-factor model versus a (k+1)-factor model is a 
 joint test whether all the free parameters (loadings) in the (k+1)st 
 column of \code{L} are zero.

 (4) For some purposes, only the subspace spanned by the factors
 is important, not the specific parameterization within this subspace.

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

  fa.ech <- echelon(fa.unrotated$loadings)
 
  fa.ech2 <- factanal(factors = 2, covmat=NetherlandsTV, rotation="echelon")
  
  cbind(loadings(fa.unrotated), loadings(fa.ech), loadings(fa.ech2))

  fa.ech3 <- echelon(fa.unrotated$loadings, reference=6:7)
  cbind(loadings(fa.unrotated), loadings(fa.ech), loadings(fa.ech3))
  }
\seealso{
   \code{\link{eiv}}, 
   \code{\link{rotations}}, 
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}}
}
\references{
   Roderick P. McDonald (1999)
     \emph{Test Theory: A Unified Treatment}, Mahwah, NJ: Erlbaum.

   Tom Wansbeek and Erik Meijer (2000) \emph{Measurement Error and 
   Latent Variables in Econometrics}, Amsterdam: North-Holland.
  }

\author{Erik Meijer and Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


