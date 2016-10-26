%%% $Id: Misra1d.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Misra1d_opt}
\alias{Misra1d}
\alias{Misra1d.f}
\alias{Misra1d.jac}
\alias{Misra1d.g}
\alias{Misra1d.setup}
\alias{Misra1d.test}
\alias{Misra1d.res}
\alias{Misra1d.h}
%%% \non_function{}
\title{Monomolecular Absorption data}
\description{
  The \code{Misra1d} data frame has 14 rows and 2 columns.
  This is the same data as \code{Misra1a} but a different model is fit.
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
plot(y ~ x, data = Misra1d)
Try(fm1 <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           start = c(b1 = 450, b2 = 0.0003) ))
Try(fm2a <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           alg = "port", start = c(b1 = 450, b2 = 0.0003) ))
Try(fm3 <- nls(y ~ b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
           start = c(b2 = 0.0005), algorithm = "plinear" ))
}
\keyword{datasets}
