\name{coef.optimx}
\alias{coef.optimx}
\encoding{UTF-8}
\title{Summarize optimx object}
\concept{minimization}
\concept{maximization}
\description{
  Summarize an \code{"optimx"} object.
}
\usage{
  \method{coef}{optimx}(object, ...) 
}
\arguments{
 \item{object}{Object returned by \code{optimx}.}
 \item{\dots}{Further arguments to be passed to the function. Currently not
   used.}
}
\value{
  \code{coef.optimx} returns the best parameters found by each method that
  returned such parameters. The returned coefficients are in the form of a
  data frame with the rows named by the relevant methods. 
}
\examples{
ans <- optimx(fn = function(x) sum(x*x), par = 1:2)
coef(ans)
}
\keyword{nonlinear}
\keyword{optimize}

