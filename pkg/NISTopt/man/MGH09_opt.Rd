%%% $Id: MGH09.Rd,v 1.4 2003/07/22 19:42:20 bates Exp $
\name{MGH09_opt}
\alias{MGH09}
\alias{MGH09.f}
\alias{MGH09.jac}
\alias{MGH09.g}
\alias{MGH09.setup}
\alias{MGH09.test}
\alias{MGH09.res}
\alias{MGH09.h}
%%% \non_function{}
\title{More, Gabrow and Hillstrom example 9}
\description{
The \code{MGH09} data frame has 11 rows and 2 columns giving
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
  This problem was found to be difficult for some very 
  good algorithms.  There is a local minimum at (+inf,
  -14.07..., -inf, -inf) with final sum of squares 
  0.00102734....
  
  See More, J. J., Garbow, B. S., and Hillstrom, K. E. 
  (1981).  \emph{Testing unconstrained optimization software.}
    \bold{ACM Transactions on Mathematical Software}. 7(1): 
    pp. 17--41.
}
\source{
Kowalik, J.S., and M. R. Osborne, (1978).  
Methods for Unconstrained Optimization Problems.  
New York, NY:  Elsevier North-Holland.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = MGH09)
## starting values for this attempt are ridiculous
Try(fm1 <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE,
           start = c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39)))
Try(fm1a <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE, alg = "port",
           start = c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39)))

Try(fm2 <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE,
           start = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39)))
Try(fm2a <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE, alg = "port",
           start = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39)))
Try(fm3 <- nls(y ~ cbind(x, x**2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE, algorithm = "plinear",
           start = c(b3 = 41.5, b4 = 39)))
Try(fm4 <- nls(y ~ cbind(x, x**2) / (x**2+x*b3+b4),
           data = MGH09, trace = TRUE, algorithm = "plinear",
           start = c(b3 = 0.415, b4 = 0.39)))
}
\keyword{datasets}
