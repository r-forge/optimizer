%%% $Id: Chwirut2.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Chwirut2_opt}
\alias{Chwirut2}
\alias{Chwirut2.f}
\alias{Chwirut2.jac}
\alias{Chwirut2.g}
\alias{Chwirut2.setup}
\alias{Chwirut2.test}
\alias{Chwirut2.res}
\alias{Chwirut2.h}


%%%\non_function{}
\title{Ultrasonic calibration data 2}
\description{
The \code{Chwirut2} data frame has nr rows and nc columns giving
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of ultrasonic response values.
    }
    \item{x}{
      A numeric vector of metal distance values.
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
  Chwirut, D., NIST (197?).  
    Ultrasonic Reference Block Study. 
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Chwirut2)
Try(fm1 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
            start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02)))
Try(fm1a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
             start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02), alg = "port"))
Try(fm2 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
            start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01)))
Try(fm2a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
             start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01), alg = "port"))
Try(fm3 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
           start = c(b1 = 0.1, p3 = 2.), alg = "plinear"))
Try(fm4 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
           start = c(b1 = 0.15, p3 = 0.01/0.008), alg = "plinear"))
}
\keyword{datasets}
