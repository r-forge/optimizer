%%% $Id: DanielWood.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{DanielWood_opt}
\alias{DanielWood}
\alias{DanielWood.f}
\alias{DanielWood.jac}
\alias{DanielWood.g}
\alias{DanielWood.setup}
\alias{DanielWood.test}
\alias{DanielWood.res}
\alias{DanielWood.h}


%%%\non_function{}
\title{Radiated energy}
\description{
The \code{DanielWood} data frame has 6 rows and 2 columns giving the
energy radiated from a carbon filament versus the absolute temperature
of the filament.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of the energy radiated from a carbon filament
      lamp.
    }
    \item{x}{
      A numeric vector of the temperature of the filament (1000 K).
    }
  }
}
\details{
    These data and model are described in Daniel and Wood
    (1980), and originally published in E.S.Keeping, 
    "Introduction to Statistical Inference," Van Nostrand
    Company, Princeton, NJ, 1962, p. 354.  The response
    variable is energy radiated from a carbon filament
    lamp per cm**2 per second, and the predictor variable
    is the absolute temperature of the filament in 1000
    degrees Kelvin.
}
\source{
Daniel, C. and F. S. Wood (1980).
Fitting Equations to Data, Second Edition. 
New York, NY:  John Wiley and Sons, pp. 428-431.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = DanielWood)
Try(fm1 <- nls(y ~ b1*x**b2, data = DanielWood, trace = TRUE,
            start = c(b1 = 1, b2 = 5)))
Try(fm1a <- nls(y ~ b1*x**b2, data = DanielWood, trace = TRUE,
            start = c(b1 = 1, b2 = 5), alg = "port"))
Try(fm2 <- nls(y ~ b1*x**b2, data = DanielWood, trace = TRUE,
            start = c(b1 = 0.7, b2 = 4)))
Try(fm2a <- nls(y ~ b1*x**b2, data = DanielWood, trace = TRUE,
            start = c(b1 = 0.7, b2 = 4), alg = "port"))
Try(fm3 <- nls(y ~ x**b2, data = DanielWood, trace = TRUE,
            start = c(b2 = 5), algorithm = "plinear"))
Try(fm4 <- nls(y ~ x**b2, data = DanielWood, trace = TRUE,
            start = c(b2 = 4), algorithm = "plinear"))
}
\keyword{datasets}
