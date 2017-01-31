\name{Test Functions}
\alias{fnTrefethen}
\alias{fnWagon}
\alias{fnColville}
\alias{fnHald}
\alias{fnAckley}
\alias{fnLangerman}
\alias{fnHartmann6}
\alias{fnMichalewicz}
\alias{fnRastrigin}
\alias{fnWhitley}
\alias{fnARWhead}
\alias{fnRana}

\title{
  Test Functions
}
\description{
  Test functions for global optimization, with several to many or even
  very many local minima and dimensions from two to ten.
}
\usage{
fnTrefethen(x)
fnWagon(x)
fnColville(x)
fnHald(x)
fnAckley(x)
fnLangerman(x)
fnHartmann6(x)
fnMichalewicz(x)
fnRastrigin(x)
fnWhitley(x)
fnARWhead(x)
fnRana(x)
}
\arguments{
  \item{x}{numeric vector, representing the real number inputs.}
}
\details{
  More details about each function can be found in the \code{testfunctions.R}
  file, like input dimension, box constraints, number of local minima, the true
  minimum and its function value, etc. This information will later on be on
  display right here.
}
\value{
  Returns a single real number, the value of a function at the given input.
}
\note{
  This package also suggests the package \emph{globalOptTests}, written by
  Kate Mullen. See below how test functions from this package can be callled.
  There are more than 50 functions, see \code{?goTest} for names.
}
\examples{
x <- rep(0.5, 5)
fnHald(x)
fnRastrigin(x)

\dontrun{
require(globalOptTests)
funBranin <- function(x) goTest(x, "Branin")
funBranin( c(0.5, 0.5) )
  }
}
