%%% $Id: Ratkowsky2.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Ratkowsky2_opt}
\alias{Ratkowsky2}
\alias{Ratkowsky2.f}
\alias{Ratkowsky2.jac}
\alias{Ratkowsky2.g}
\alias{Ratkowsky2.setup}
\alias{Ratkowsky2.test}
\alias{Ratkowsky2.res}
\alias{Ratkowsky2.h}
%%% \non_function{}
\title{Pasture yield data}
\description{
  The \code{Ratkowsky2} data frame has 9 rows and 2 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of pasture yields.
    }
    \item{x}{
      A numeric vector of growing times.
    }
  }
}
\details{
  This model and data are an example of fitting
  sigmoidal growth curves taken from Ratkowsky (1983).
  The response variable is pasture yield, and the
  predictor variable is growing time.
}
\source{
Ratkowsky, D.A. (1983).  
Nonlinear Regression Modeling.
New York, NY:  Marcel Dekker, pp. 61 and 88.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Ratkowsky2)

Try(fm1 <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
           start = c(b1 = 100, b2 = 1, b3 = 0.1)))
Try(fm1a <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2,
            trace = TRUE, alg = "port", 
            start = c(b1 = 100, b2 = 1, b3 = 0.1)))
Try(fm2 <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
           start = c(b1 = 75, b2 = 2.5, b3 = 0.07)))
Try(fm2a <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2,
            trace = TRUE, alg = "port", 
            start = c(b1 = 75, b2 = 2.5, b3 = 0.07)))
Try(fm3 <- nls(y ~ 1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
           start = c(b2 = 1, b3 = 0.1), alg = "plinear"))
Try(fm4 <- nls(y ~ 1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
           start = c(b2 = 2.5, b3 = 0.07), alg = "plinear"))

## Using a self-starting model
Try(fm5 <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = Ratkowsky2))
summary(fm5)
}         
\keyword{datasets}
