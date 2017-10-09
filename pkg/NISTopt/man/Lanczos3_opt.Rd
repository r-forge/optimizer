%%% $Id: Lanczos3.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Lanczos3_opt}
\alias{Lanczos3}
\alias{Lanczos3.f}
\alias{Lanczos3.jac}
\alias{Lanczos3.g}
\alias{Lanczos3.setup}
\alias{Lanczos3.test}
\alias{Lanczos3.res}
\alias{Lanczos3.h}
%%% \non_function{}
\title{Generated data}
\description{
The \code{Lanczos3} data frame has 24 rows and 2 columns of generated data.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of generated responses.
    }
    \item{x}{
      A numeric vector of generated input values.
    }
  }
}
\details{
  These data are taken from an example discussed in
  Lanczos (1956).  The data were generated to 5-digits
  of accuracy using
  \code{f(x) = 0.0951*exp(-x) + 0.8607*exp(-3*x) + 1.5576*exp(-5*x)}.
}
\source{
    Lanczos, C. (1956).
    Applied Analysis.
    Englewood Cliffs, NJ:  Prentice Hall, pp. 272-280.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Lanczos3)
## plot log response to see the number of exponential terms
plot(y ~ x, data = Lanczos3, log = "y")
Try(fm1 <- nls(y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x),
           data = Lanczos3, trace = TRUE,
           start = c(b1 = 1.2, b2 = 0.3, b3 = 5.6, b4 = 5.5,
                     b5 = 6.5, b6 = 7.6)))
Try(fm1a <- nls(y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x),
           data = Lanczos3, trace = TRUE, alg = "port",
           start = c(b1 = 1.2, b2 = 0.3, b3 = 5.6, b4 = 5.5,
                     b5 = 6.5, b6 = 7.6)))
Try(fm2 <- nls(y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x),
           data = Lanczos3, trace = TRUE,
           start = c(b1 = 0.5, b2 = 0.7, b3 = 3.6, b4 = 4.2,
                     b5 = 4, b6 = 6.3)))
Try(fm2a <- nls(y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x),
           data = Lanczos3, trace = TRUE, alg = "port",
           start = c(b1 = 0.5, b2 = 0.7, b3 = 3.6, b4 = 4.2,
                     b5 = 4, b6 = 6.3)))
Try(fm3 <- nls(y ~ exp(outer(x,-c(b2, b4, b6))),
           data = Lanczos3, trace = TRUE, algorithm = "plinear",
           start = c(b2 = 0.3, b4 = 5.5, b6 = 7.6)))
Try(fm4 <- nls(y ~ exp(outer(x,-c(b2, b4, b6))),
           data = Lanczos3, trace = TRUE, algorithm = "plinear",
           start = c(b2 = 0.7, b4 = 4.2, b6 = 6.3)))
## Use analytic derivatives
Lanczos <- deriv(~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x),
                 paste("b", 1:6, sep = ""),
                 function(x, b1, b2, b3, b4, b5, b6){})
Try(fm5 <- nls(y ~ Lanczos(x, b1, b2, b3, b4, b5, b6),
           data = Lanczos3, trace = TRUE,
           start = c(b1 = 1.2, b2 = 0.3, b3 = 5.6, b4 = 5.5,
                     b5 = 6.5, b6 = 7.6)))
Try(fm5a <- nls(y ~ Lanczos(x, b1, b2, b3, b4, b5, b6),
           data = Lanczos3, trace = TRUE, alg = "port",
           start = c(b1 = 1.2, b2 = 0.3, b3 = 5.6, b4 = 5.5,
                     b5 = 6.5, b6 = 7.6)))
Try(fm6 <- nls(y ~ Lanczos(x, b1, b2, b3, b4, b5, b6),
           data = Lanczos3, trace = TRUE,
           start = c(b1 = 0.5, b2 = 0.7, b3 = 3.6, b4 = 4.2,
                     b5 = 4, b6 = 6.3)))
Try(fm6a <- nls(y ~ Lanczos(x, b1, b2, b3, b4, b5, b6),
           data = Lanczos3, trace = TRUE, alg = "port",
           start = c(b1 = 0.5, b2 = 0.7, b3 = 3.6, b4 = 4.2,
                     b5 = 4, b6 = 6.3)))
}
\keyword{datasets}
