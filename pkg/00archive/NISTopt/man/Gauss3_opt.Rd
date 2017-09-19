%%% $Id: Gauss3.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Gauss3_opt}
\alias{Gauss3}
\alias{Gauss3.f}
\alias{Gauss3.jac}
\alias{Gauss3.g}
\alias{Gauss3.setup}
\alias{Gauss3.test}
\alias{Gauss3.res}
\alias{Gauss3.h}
%%% \non_function{}
\title{Generated data}
\description{
The \code{Gauss3} data frame has 250 rows and 2 columns giving generated
data of Gaussian peaks with a decaying exponential background.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of generated responses.
    }
    \item{x}{
      A numeric vector of generated inputs.
    }
  }
}
\details{
  The data are two strongly-blended Gaussians on a 
  decaying exponential baseline plus normally 
  distributed zero-mean noise with variance = 6.25.
}
\source{
Rust, B., NIST (1996).
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Gauss3)
Try(fm1 <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss3, trace = TRUE,
           start = c(b1 = 94.9, b2 = 0.009, b3 = 90.1, b4 = 113, b5 = 20,
                     b6 = 73.8, b7 = 140, b8 = 20)))
Try(fm1a <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
                + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss3, trace = TRUE,
            start = c(b1 = 94.9, b2 = 0.009, b3 = 90.1, b4 = 113, b5 = 20,
                      b6 = 73.8, b7 = 140, b8 = 20), alg = "port"))
Try(fm2 <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss3, trace = TRUE,
           start = c(b1 = 96, b2 = 0.0096, b3 = 80, b4 = 110, b5 = 25,
                     b6 = 74, b7 = 139, b8 = 25)))
Try(fm2a <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss3, trace = TRUE,
           start = c(b1 = 96, b2 = 0.0096, b3 = 80, b4 = 110, b5 = 25,
                     b6 = 74, b7 = 139, b8 = 25), alg = "port"))
Try(fm3 <- nls(y ~ cbind(exp(-b2*x), exp(-(x-b4)**2/b5**2), exp(-(x-b7)**2/b8**2)),
           data = Gauss3, trace = TRUE,
           start = c(b2 = 0.009, b4 = 113, b5 = 20, b7 = 140, b8 = 20),
           algorithm = "plinear"))
Try(fm4 <- nls(y ~ cbind(exp(-b2*x), exp(-(x-b4)**2/b5**2), exp(-(x-b7)**2/b8**2)),
           data = Gauss3, trace = TRUE,
           start = c(b2 = 0.0096, b4 = 110, b5 = 25, b7 = 139, b8 = 25),
           algorithm = "plinear"))
}
\keyword{datasets}
