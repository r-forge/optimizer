%%% $Id: Eckerle4.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Eckerle4_opt}
\alias{Eckerle4}
\alias{Eckerle4.f}
\alias{Eckerle4.jac}
\alias{Eckerle4.g}
\alias{Eckerle4.setup}
\alias{Eckerle4.test}
\alias{Eckerle4.res}
\alias{Eckerle4.h}


%%%\non_function{}
\title{Circular interference data}
\description{
The \code{Eckerle4} data frame has 35 rows and 2 columns giving
transmittance as a function of wavelength.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of transmittance values.
    }
    \item{x}{
      A numeric vector of wavelengths.
    }
  }
}
\details{
  These data are the result of a NIST study involving
  circular interference transmittance.  The response
  variable is transmittance, and the predictor variable
  is wavelength.
}
\source{
Eckerle, K., NIST (197?).  
Circular Interference Transmittance Study.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Eckerle4)
## should fail - ridiculous starting value for b3
Try(fm1 <- nls(y ~ (b1/b2) * exp(-0.5*((x-b3)/b2)**2), Eckerle4,
               trace = TRUE,
               start = c(b1 = 1, b2 = 10, b3 = 500)))
Try(fm1a <- nls(y ~ (b1/b2) * exp(-0.5*((x-b3)/b2)**2), Eckerle4,
                  trace = TRUE, alg = "port",
                  start = c(b1 = 1, b2 = 10, b3 = 500)))
Try(fm2 <- nls(y ~ (b1/b2) * exp(-0.5*((x-b3)/b2)**2),
            Eckerle4, trace = TRUE,
            start = c(b1 = 1.5, b2 = 5, b3 = 450)))
Try(fm2a <- nls(y ~ (b1/b2) * exp(-0.5*((x-b3)/b2)**2),
            Eckerle4, trace = TRUE, alg = "port",
            start = c(b1 = 1.5, b2 = 5, b3 = 450)))
## should fail - ridiculous starting value for b3
Try(fm3 <- nls(y ~ (1/b2) * exp(-0.5*((x-b3)/b2)**2),
               Eckerle4, trace = TRUE,
               start = c(b2 = 10, b3 = 500), algorithm = "plinear"))
Try(fm4 <- nls(y ~ (1/b2) * exp(-0.5*((x-b3)/b2)**2), Eckerle4, trace = TRUE,
           start = c(b2 = 5, b3 = 450), algorithm = "plinear"))
}
\keyword{datasets}
