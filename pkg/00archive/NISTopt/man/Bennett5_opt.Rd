%%% $Id: Bennett5.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Bennett5_opt}
\alias{Bennett5}
\alias{Bennett5.f}
\alias{Bennett5.jac}
\alias{Bennett5.g}
\alias{Bennett5.setup}
\alias{Bennett5.test}
\alias{Bennett5.res}
\alias{Bennett5.h}
%%% \non_function{}
\title{Magentization modelling}
\description{
  Data and functions for fitting the Bennett5 nonlinear least squares problem.
  The \code{Bennett5} data frame has 154 rows and 2 columns of data from a
  magnetism study
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of magnetism values.
    }
    \item{x}{
      A numeric vector of log(time).
    }
  }
}
\usage{
Bennett5.f(b)
}
\arguments{
 \item{b}{A numeric vector of starting estimates.}
}
\value{
 The numeric value of the sum of squared deviations between model and data.
}
\details{
  These data are the result of a NIST study involving superconductivity
  magnetization modeling.  The response variable is magnetism, and the
  predictor variable is the log of time in minutes.
}
\source{
  Bennett, L., L. Swartzendruber, and H. Brown, 
  NIST (1994).  
  Superconductivity Magnetization Modeling.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = NISTopt:::Bennett5)
Try(fm1 <- nls(y ~ b1*(b2+x)**(-1/b3), data = NISTopt:::Bennett5,
           start = c(b1 = -2000, b2 = 50, b3 = 0.8), trace = TRUE))
Try(fm1a <- nls(y ~ b1*(b2+x)**(-1/b3), data = NISTopt:::Bennett5,
            start = c(b1 = -2000, b2 = 50, b3 = 0.8),
            trace = TRUE, alg = "port"))
Try(fm2 <- nls(y ~ b1*(b2+x)**(-1/b3), data = NISTopt:::Bennett5,
           start = c(b1 = -1500, b2 = 45, b3 = 0.85), trace = TRUE))
Try(fm2a <- nls(y ~ b1*(b2+x)**(-1/b3), data = NISTopt:::Bennett5,
            start = c(b1 = -1500, b2 = 45, b3 = 0.85),
            trace = TRUE, alg = "port"))
Try(fm3 <- nls(y ~ (b2+x)**(-1/b3), data = NISTopt:::Bennett5, alg = "plinear",
           start = c( b2 = 50, b3 = 0.8), trace = TRUE))
Try(fm4 <- nls(y ~ (b2+x)**(-1/b3), data = NISTopt:::Bennett5, alg = "plinear",
           start = c( b2 = 45, b3 = 0.8), trace = TRUE))
}
\keyword{datasets}
