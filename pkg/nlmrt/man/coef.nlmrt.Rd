\name{coef.nlmrt}
\alias{coef.nlmrt}
\title{Output model coefficients for nlmrt object.
}
\description{
    \code{coef.nlmrt} extracts and displays the coefficients for a model 
    estimated by \code{nlxb} or \code{nlfb} in the \code{nlmrt} structured
    \code{object}. 
}
\usage{
    \method{coef}{nlmrt}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
     An object of class 'nlmrt'
}
  \item{\dots}{
     Any data needed for the function. We do not know of any!
}
}
\details{
    \code{coef.nlmrt} extracts and displays the coefficients for a model 
    estimated by \code{nlxb} or \code{nlfb}. 
}
\value{
    returns the coefficients from the nlmrt object.
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
