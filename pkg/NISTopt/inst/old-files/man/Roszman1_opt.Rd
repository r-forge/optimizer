%%% $Id: Roszman1.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Roszman1_opt}
\alias{Roszman1}
\alias{Roszman1.f}
\alias{Roszman1.jac}
\alias{Roszman1.g}
\alias{Roszman1.setup}
\alias{Roszman1.test}
\alias{Roszman1.res}
\alias{Roszman1.h}
%%% \non_function{}
\title{Quantum defects in iodine}
\description{
The \code{Roszman1} data frame has 25 rows and 2 columns of data on the
number of quantum defects in iodine atoms at different energy states.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of number of quantum defects.
    }
    \item{x}{
      A numeric vector of the excited energy state.
    }
  }
}
\details{
These data are the result of a NIST study involving
quantum defects in iodine atoms.  The response
variable is the number of quantum defects, and the
predictor variable is the excited energy state.
The argument to the ARCTAN function is in radians.
}
\source{
Roszman, L., NIST (19??).  
Quantum Defects for Sulfur I Atom.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Roszman1)
Try(fm1 <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
           start = c(b1 = 0.1, b2 = -0.00001, b3 = 1000, b4 = -100),
           trace = TRUE))
Try(fm1a <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
           start = c(b1 = 0.1, b2 = -0.00001, b3 = 1000, b4 = -100),
           alg = "port", trace = TRUE))
Try(fm2 <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
           start = c(b1 = 0.2, b2 = -0.0000015, b3 = 1200, b4 = -150),
           trace = TRUE))
Try(fm2a <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
           start = c(b1 = 0.2, b2 = -0.0000015, b3 = 1200, b4 = -150),
           alg = "port", trace = TRUE))
}
\keyword{datasets}
