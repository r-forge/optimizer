%%% $Id: Misra1b.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Misra1b_opt}
\alias{Misra1b}
\alias{Misra1b.f}
\alias{Misra1b.jac}
\alias{Misra1b.g}
\alias{Misra1b.setup}
\alias{Misra1b.test}
\alias{Misra1b.res}
\alias{Misra1b.h}
%%% \non_function{}
\title{Monomolecular Absorption Data}
\description{
The \code{Misra1b} data frame has 14 rows and 2 columns.
It is the same data as \code{Misra1a} but a different model is fit.
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
plot(y ~ x, data = Misra1b)
Try(fm1 <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
           start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
           alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
           start = c(b1 = 300, b2 = 0.0002) ))
Try(fm2a <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
           alg = "port", start = c(b1 = 300, b2 = 0.0002) ))
Try(fm3 <- nls(y ~ 1-(1+b2*x/2)**(-2), data = Misra1b, trace = TRUE,
           start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ 1-(1+b2*x/2)**(-2), data = Misra1b, trace = TRUE,
           start = c(b2 = 0.0005), algorithm = "plinear" ))
}
\keyword{datasets}
