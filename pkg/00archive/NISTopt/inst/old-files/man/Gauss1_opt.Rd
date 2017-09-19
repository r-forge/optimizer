%%% $Id: Gauss1.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Gauss1_opt}
\alias{Gauss1}
\alias{Gauss1.f}
\alias{Gauss1.jac}
\alias{Gauss1.g}
\alias{Gauss1.setup}
\alias{Gauss1.test}
\alias{Gauss1.res}
\alias{Gauss1.h}
%%% \non_function{}
\title{Generated data}
\description{
The \code{Gauss1} data frame has 250 rows and 2 columns of generated data.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of generated responses.
    }
    \item{x}{
      A numeric vector of generated input values.
    }
  }
}
\details{
    The data are generated data with two well-separated Gaussians on a 
    decaying exponential baseline plus normally 
    distributed zero-mean noise with variance = 6.25.
}
\source{
    Rust, B., NIST (1996).
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Gauss1)
Try(fm1 <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss1, trace = TRUE,
            start = c(b1 = 97.0, b2 = 0.009, b3 = 100.0, b4 = 65.0, b5 = 20.0,
                     b6 = 70.0, b7 = 178., b8 = 16.5)))
Try(fm1a <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
                + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss1, trace = TRUE,
             start = c(b1 = 97.0, b2 = 0.009, b3 = 100.0, b4 = 65.0, b5 = 20.0,
                     b6 = 70.0, b7 = 178., b8 = 16.5), alg = "port"))
Try(fm2 <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss1, trace = TRUE,
           start = c(b1 = 94.0, b2 =  0.0105, b3 = 99.0, b4 = 63.0, b5 = 25.0,
                     b6 = 71.0, b7 = 180.0, b8 = 20.0)))
Try(fm2a <- nls(y ~ b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
               + b6*exp( -(x-b7)**2 / b8**2 ), data = Gauss1, trace = TRUE,
           start = c(b1 = 94.0, b2 =  0.0105, b3 = 99.0, b4 = 63.0, b5 = 25.0,
                     b6 = 71.0, b7 = 180.0, b8 = 20.0), alg = "port"))
Try(fm3 <- nls(y ~ cbind(exp(-b2*x), exp(-(x-b4)**2/b5**2), exp(-(x-b7)**2/b8**2)),
           data = Gauss1, trace = TRUE,
           start = c( b2 = 0.009, b4 = 65.0, b5 = 20.0, b7 = 178., b8 = 16.5),
           algorithm = "plinear"))
Try(fm4 <- nls(y ~ cbind(exp(-b2*x), exp(-(x-b4)**2/b5**2), exp(-(x-b7)**2/b8**2)),
           data = Gauss1, trace = TRUE,
           start = c( b2 = 0.0105, b4 = 63.0, b5 = 25.0, b7 = 180., b8 = 20.0),
           algorithm = "plinear"))
}
\keyword{datasets}
