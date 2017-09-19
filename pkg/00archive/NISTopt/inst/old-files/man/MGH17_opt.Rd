%%% $Id: MGH17.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{MGH17_opt}
\alias{MGH17}
\alias{MGH17.f}
\alias{MGH17.jac}
\alias{MGH17.g}
\alias{MGH17.setup}
\alias{MGH17.test}
\alias{MGH17.res}
\alias{MGH17.h}
%%% \non_function{}
\title{More, Gabrow and Hillstrom example 17}
\description{
The \code{MGH17} data frame has 33 rows and 2 columns
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of response values.
    }
    \item{x}{
      A numeric vector of input values.
    }
  }
}
\details{
  This problem was found to be difficult for some very
  good algorithms.

  See More, J. J., Garbow, B. S., and Hillstrom, K. E.
  (1981).  \emph{Testing unconstrained optimization software.}
  \bold{ACM Transactions on Mathematical Software.} 7(1):
  pp. 17-41.
}
\source{
  Osborne, M. R. (1972).  
  Some aspects of nonlinear least squares 
  calculations.  In Numerical Methods for Nonlinear 
  Optimization, Lootsma (Ed).  
  New York, NY:  Academic Press, pp. 171-189.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = MGH17)

## Starting values here are ridiculous
Try(fm1 <- nls(y ~ b1 + b2*exp(-x*b4) + b3*exp(-x*b5),
           data = MGH17, trace = TRUE,
           start = c(b1 = 50, b2 = 150, b3 = -100, b4 = 1, b5 = 2)))
Try(fm1a <- nls(y ~ b1 + b2*exp(-x*b4) + b3*exp(-x*b5),
           data = MGH17, trace = TRUE, alg = "port",
           start = c(b1 = 50, b2 = 150, b3 = -100, b4 = 1, b5 = 2)))

Try(fm2 <- nls(y ~ b1 + b2*exp(-x*b4) + b3*exp(-x*b5),
           data = MGH17, trace = TRUE,
           start = c(b1 = 0.5, b2 = 1.5, b3 = -1, b4 = 0.01, b5 = 0.02)))
Try(fm2a <- nls(y ~ b1 + b2*exp(-x*b4) + b3*exp(-x*b5),
           data = MGH17, trace = TRUE, alg = "port",
           start = c(b1 = 0.5, b2 = 1.5, b3 = -1, b4 = 0.01, b5 = 0.02)))

Try(fm3 <- nls(y ~ cbind(1, exp(-x*b4), exp(-x*b5)),
           data = MGH17, trace = TRUE, algorithm = "plinear",
           start = c(b4 = 1, b5 = 2)))

Try(fm4 <- nls(y ~ cbind(1, exp(-x*b4), exp(-x*b5)),
           data = MGH17, trace = TRUE, algorithm = "plinear",
           start = c(b4 = 0.01, b5 = 0.02)))
}
\keyword{datasets}
