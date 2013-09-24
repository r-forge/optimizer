\name{print.nlmrt}
\alias{print.nlmrt}
\title{Print method for an object of class nlmrt.
}
\description{
    Print summary output (but involving some serious computations!) of
    an object of class nlmrt from \code{nlxb} or \code{nlfb} from package
    \code{nlmrt}.
}
\usage{
    \method{print}{nlmrt}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An object of class 'nlmrt'
}
  \item{\dots}{
     Any data needed for the function. We do not know of any!
}
}
\details{
    \code{printsum.nlmrt} performs a print method for an object of class 'nlmrt' that
    has been created by a routine such as \code{nlfb} or \code{nlxb} for nonlinear
    least squares problems.
}
\value{
    Invisibly returns the input object.
}
\author{
John C Nash <nashjc@uottawa.ca>
}
\note{
   Special notes, if any, will appear here.
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\keyword{ nonlinear least squares }

