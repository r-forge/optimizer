%%% $Id: Chwirut1.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Chwirut1_opt}
\alias{Chwirut1}
\alias{Chwirut1.f}
\alias{Chwirut1.jac}
\alias{Chwirut1.g}
\alias{Chwirut1.setup}
\alias{Chwirut1.test}
\alias{Chwirut1.res}
\alias{Chwirut1.h}

%%%\non_function{}
\title{Ultrasonic calibration study 1}
\description{
  The \code{Chwirut1} data frame has 214 rows and 2 columns giving
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of ultrasonic response values
    }
    \item{x}{
      A numeric vector or metal distance values
    }
  }
}
\details{
  These data are the result of a NIST study involving
  ultrasonic calibration.  The response variable is
  ultrasonic response, and the predictor variable is
  metal distance.
}
\source{
  Chwirut, D., NIST (197?).   Ultrasonic Reference Block Study. 
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Chwirut1)
Try(fm1 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
           start = c(b1 = 0.1, b2 = 0.01, b3 = 0.02)))
Try(fm1a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
           start = c(b1 = 0.1, b2 = 0.01, b3 = 0.02), alg = "port"))
Try(fm2 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
           start = c(b1 = 0.15, b2 = 0.008, b3 = 0.010)))
Try(fm2a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
            start = c(b1 = 0.15, b2 = 0.008, b3 = 0.010), alg = "port"))
Try(fm3 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut1, trace = TRUE,
           start = c(b1 = 0.1, p3 = 0.02/0.01), algorithm = "plinear"))
Try(fm4 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut1, trace = TRUE,
           start = c(b1 = 0.15, p3 = 0.01/0.008), algorithm = "plinear"))
}
\keyword{datasets}
