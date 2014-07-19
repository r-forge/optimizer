\name{print.nls14}
\alias{print.nls14}
\title{Print method for an object of class nls14.
}
\description{
    Print summary output (but involving some serious computations!) of
    an object of class nls14 from \code{nlxb} or \code{nlfb} from package
    \code{nls14}.
}
\usage{
    \method{print}{nls14}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An object of class 'nls14'
}
  \item{\dots}{
     Any data needed for the function. We do not know of any!
}
}
\details{
    \code{printsum.nls14} performs a print method for an object of class 'nls14' that
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

