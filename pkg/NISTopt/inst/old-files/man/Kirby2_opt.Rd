%%% $Id: Kirby2.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{Kirby2_opt}
\alias{Kirby2}
\alias{Kirby2.f}
\alias{Kirby2.jac}
\alias{Kirby2.g}
\alias{Kirby2.setup}
\alias{Kirby2.test}
\alias{Kirby2.res}
\alias{Kirby2.h}
%%% \non_function{}
\title{Microscope line width standards}
\description{
The \code{Kirby2} data frame has 151 rows and 2 columns of data from an
NIST study on scanning electron microscope line width standards.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of response values.
    }
    \item{x}{
      A numeric vector of input values.
    }
  }
}
\details{
These data are the result of a NIST study involving
scanning electron microscope line with standards.
}
\source{
Kirby, R., NIST (197?).  
Scanning electron microscope line width standards.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Kirby2)
Try(fm1 <- nls(y ~ (b1 + b2*x + b3*x**2) / (1 + b4*x + b5*x**2),
           data = Kirby2, trace = TRUE,
           start = c(b1 = 2, b2 = -0.1, b3 = 0.003,
                     b4 = -0.001, b5 = 0.00001)))
Try(fm1a <- nls(y ~ (b1 + b2*x + b3*x**2) / (1 + b4*x + b5*x**2),
           data = Kirby2, trace = TRUE, alg = "port",
           start = c(b1 = 2, b2 = -0.1, b3 = 0.003,
                     b4 = -0.001, b5 = 0.00001)))
Try(fm2 <- nls(y ~ (b1 + b2*x + b3*x**2) / (1 + b4*x + b5*x**2),
           data = Kirby2, trace = TRUE,
           start = c(b1 = 1.5, b2 = -0.15, b3 = 0.0025,
                     b4 = -0.0015, b5 = 0.00002)))
Try(fm2a <- nls(y ~ (b1 + b2*x + b3*x**2) / (1 + b4*x + b5*x**2),
             data = Kirby2, trace = TRUE, alg = "port",
             start = c(b1 = 1.5, b2 = -0.15, b3 = 0.0025,
                       b4 = -0.0015, b5 = 0.00002)))
Try(fm3 <- nls(y ~ cbind(1, x, x**2)/(1 + x*(b4 + b5*x)),
           data = Kirby2, trace = TRUE, algorithm = "plinear",
           start = c(b4 = -0.001, b5 = 0.00001)))
Try(fm4 <- nls(y ~ cbind(1, x, x**2)/(1 + x*(b4 + b5*x)),
           data = Kirby2, trace = TRUE, algorithm = "plinear",
           start = c(b4 = -0.0015, b5 = 0.00002)))
}
\keyword{datasets}
