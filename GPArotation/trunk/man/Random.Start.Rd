\name{Random.Start}
\alias{Random.Start}

\title{Generate a Random Orthogonal Rotation}
\usage{
    Random.Start(k)
    }
\arguments{
    \item{k}{An integer indicating the dimension of the square matrix.}
    }

\description{
   Random orthogonal rotation to use as Tmat matrix to start GPForth or GPFoblq.
   }
\value{An orthogonal matrix.}

\details{
  The random start function produces an orthogonal matrix with columns
  of length one based on the QR decompostion.
  }
\seealso{
   \code{\link{GPForth}}
   \code{\link{GPFoblq}}
   \code{\link{oblimin}}
   }
\examples{
   Global.min <- function(A,method,B=10){
      fv <- rep(0,B)
      seeds <- sample(1e+7, B)
      for(i in 1:B){
    	cat(i," ")
    	set.seed(seeds[i])
    	gpout <- GPFoblq(A=A, Random.Start(ncol(A)), method=method)
    	dtab <- dim(gpout$Table)
    	fv[i] <- gpout$Table[dtab[1],2]
    	cat(fv[i], "\n")
      }
      cat("Min is ",min(fv),"\n")
      set.seed(seeds[order(fv)[1]])
      ans <- GPFoblq(A=A, Random.Start(ncol(A)), method=method)
      ans
      }

   data("Thurstone", package="GPArotation")

   Global.min(box26,"simplimax",10)
   }
\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert
	}
\concept{rotation}
\keyword{multivariate}


