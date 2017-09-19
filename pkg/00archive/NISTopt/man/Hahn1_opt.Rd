%%% $Id: Hahn1.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Hahn1_opt}
\alias{Hahn1}
\alias{Hahn1.f}
\alias{Hahn1.jac}
\alias{Hahn1.g}
\alias{Hahn1.setup}
\alias{Hahn1.test}
\alias{Hahn1.res}
\alias{Hahn1.h}
%%% \non_function{}
\title{Thermal expansion data}
\description{
The \code{Hahn1} data frame has 236 rows and 2 columns of data from a
study on the thermal expansion of copper.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of values of the coefficient of thermal expansion.
    }
    \item{x}{
      A numeric vector of temperatures (K).
    }
  }
}
\details{
  These data are the result of a NIST study involving
  the thermal expansion of copper.  The response 
  variable is the coefficient of thermal expansion, and
  the predictor variable is temperature in degrees 
  Kelvin.
}
\source{
Hahn, T., NIST (197?). 
Copper Thermal Expansion Study.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Hahn1)
Try(fm1 <- nls(y ~ (b1+b2*x+b3*x**2+b4*x**3) / (1+b5*x+b6*x**2+b7*x**3),
           data = Hahn1, trace = TRUE,
           start = c(b1 = 10, b2 = -1, b3 = 0.05,
             b4 = -0.00001, b5 = -0.05, b6 = 0.001, b7 = -0.000001)))
Try(fm1a <- nls(y ~ (b1+b2*x+b3*x**2+b4*x**3) / (1+b5*x+b6*x**2+b7*x**3),
           data = Hahn1, trace = TRUE, alg = "port",
           start = c(b1 = 10, b2 = -1, b3 = 0.05,
             b4 = -0.00001, b5 = -0.05, b6 = 0.001, b7 = -0.000001)))
Try(fm2 <- nls(y ~ (b1+b2*x+b3*x**2+b4*x**3) / (1+b5*x+b6*x**2+b7*x**3),
           data = Hahn1, trace = TRUE,
           start = c(b1 = 1, b2 = -0.1, b3 = 0.005, b4 = -0.000001,
              b5 = -0.005, b6 = 0.0001, b7 = -0.0000001)))
Try(fm2a <- nls(y ~ (b1+b2*x+b3*x**2+b4*x**3) / (1+b5*x+b6*x**2+b7*x**3),
             data = Hahn1, trace = TRUE, alg = "port",
           start = c(b1 = 1, b2 = -0.1, b3 = 0.005, b4 = -0.000001,
              b5 = -0.005, b6 = 0.0001, b7 = -0.0000001)))
Try(fm3 <- nls(y ~ cbind(1, x, x^2, x^3)/(1+x*(b5+x*(b6+x*b7))),
           data = Hahn1, trace = TRUE, algorithm = "plinear",
           start = c(b5 = -0.05, b6 = 0.001, b7 = -0.000001)))
Try(fm4 <- nls(y ~ cbind(1, x, x^2, x^3)/(1+x*(b5+x*(b6+x*b7))),
           data = Hahn1, trace = TRUE, algorithm = "plinear",
           start = c(b5 = -0.005, b6 = 0.0001, b7 = -0.0000001)))
}
\keyword{datasets}
