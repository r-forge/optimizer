%%% $Id: Ratkowsky3.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{Ratkowsky3_opt}
\alias{Ratkowsky3}
\alias{Ratkowsky3.f}
\alias{Ratkowsky3.jac}
\alias{Ratkowsky3.g}
\alias{Ratkowsky3.setup}
\alias{Ratkowsky3.test}
\alias{Ratkowsky3.res}
\alias{Ratkowsky3.h}
%%% \non_function{}
\title{Onion growth data}
\description{
The \code{Ratkowsky3} data frame has 15 rows and 2 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of dry weights of onion bulbs and tops.
    }
    \item{x}{
      A numeric vector of growing times.
    }
  }
}
\details{
This model and data are an example of fitting  
sigmoidal growth curves taken from Ratkowsky (1983).  
The response variable is the dry weight of onion bulbs 
and tops, and the predictor variable is growing time. 
}
\source{
Ratkowsky, D.A. (1983).  
Nonlinear Regression Modeling.
New York, NY:  Marcel Dekker, pp. 62 and 88.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Ratkowsky3)

## causes NA/NaN/Inf error
Try(fm1 <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b1 = 100, b2 = 10, b3 = 1, b4 = 1),
           trace = TRUE))
Try(fm1a <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b1 = 100, b2 = 10, b3 = 1, b4 = 1),
           alg = "port", trace = TRUE))

Try(fm2 <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b1 = 700, b2 = 5, b3 = 0.75, b4 = 1.3),
           trace = TRUE))
Try(fm2a <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b1 = 700, b2 = 5, b3 = 0.75, b4 = 1.3),
           alg = "port", trace = TRUE))

Try(fm3 <- nls(y ~ 1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b2 = 10, b3 = 1, b4 = 1), algorithm = "plinear",
           trace = TRUE))
Try(fm4 <- nls(y ~ 1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
           start = c(b2 = 5, b3 = 0.75, b4 = 1.3), algorithm = "plinear",
           trace = TRUE))
}
\keyword{datasets}
