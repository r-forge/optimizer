%%% $Id: Nelson.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Nelson_opt}
\alias{Nelson}
\alias{Nelson.f}
\alias{Nelson.jac}
\alias{Nelson.g}
\alias{Nelson.setup}
\alias{Nelson.test}
\alias{Nelson.res}
\alias{Nelson.h}
%%% \non_function{}
\title{Dialectric breakdown data}
\description{
The \code{Nelson} data frame has 128 rows and 3 columns of data from an
accelerated test of dialectric breakdown.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of dialectric breakdown strength values.
    }
    \item{x1}{
      A numeric vector of time values.
    }
    \item{x2}{
      A numeric vector of temperature values.
    }
  }
}
\details{
  These data are the result of a study involving
  the analysis of performance degradation data from
  accelerated tests, published in IEEE Transactions
  on Reliability.  The response variable is dialectric
  breakdown strength in kilo-volts, and the predictor
  variables are time in weeks and temperature in degrees
  Celsius.
}
\source{
Nelson, W. (1981).  
Analysis of Performance-Degradation Data.  
IEEE Transactions on Reliability.
Vol. 2, R-30, No. 2, pp. 149-155.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x1, data = Nelson, log = "y")
plot(y ~ x2, data = Nelson, log = "y")
coplot(y ~ x1 | x2, data = Nelson)
coplot(y ~ x2 | x1, data = Nelson)

Try(fm1 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
           start = c(b1 = 2, b2 = 0.0001, b3 = -0.01), trace = TRUE))
Try(fm1a <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
            trace = TRUE, alg = "port",
            start = c(b1 = 2, b2 = 0.0001, b3 = -0.01)))

Try(fm2 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
           start = c(b1 = 2.5, b2 = 0.000000005, b3 = -0.05), trace = TRUE))
Try(fm2 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
           trace = TRUE, alg = "port", 
           start = c(b1 = 2.5, b2 = 0.000000005, b3 = -0.05)))

Try(fm3 <- nls(log(y) ~ cbind(1, -x1 * exp(-b3*x2)), data = Nelson,
           start = c(b3 = -0.01), trace = TRUE, algorithm = "plinear"))

Try(fm4 <- nls(log(y) ~ cbind(1, -x1 * exp(-b3*x2)), data = Nelson,
           start = c(b3 = -0.05), trace = TRUE, algorithm = "plinear"))
}
\keyword{datasets}
