%%% $Id: Misra1a.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Misra1a_opt}
\alias{Misra1a}
\alias{Misra1a.f}
\alias{Misra1a.jac}
\alias{Misra1a.g}
\alias{Misra1a.setup}
\alias{Misra1a.test}
\alias{Misra1a.res}
\alias{Misra1a.h}
%%% \non_function{}
\title{Monomolecular Absorption Data}
\description{
  The \code{Misra1a} data frame has 14 rows and 2 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of volume values.
    }
    \item{x}{
      A numeric vector of pressure values.
    }
  }
}
\details{
  These data are the result of a NIST study regarding
  dental research in monomolecular adsorption.  The
  response variable is volume, and the predictor
  variable is pressure.
}
\source{
Misra, D., NIST (1978).
Dental Research Monomolecular Adsorption Study.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Misra1a)
Try(fm1 <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
           start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
           alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
           start = c(b1 = 250, b2 = 0.0005) ))
Try(fm2a <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
           alg = "port", start = c(b1 = 250, b2 = 0.0005) ))
Try(fm3 <- nls(y ~ 1-exp(-b2*x), data = Misra1a, trace = TRUE,
           start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ 1-exp(-b2*x), data = Misra1a, trace = TRUE,
           start = c(b2 = 0.0005), algorithm = "plinear" ))

## Using a self-starting model
Try(fm5 <- nls(y ~ SSasympOrig(x, Asym, lrc), data = Misra1a))
}
\keyword{datasets}
