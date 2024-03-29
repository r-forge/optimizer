% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.nlsr.R
\name{coef.nlsr}
\alias{coef.nlsr}
\title{coef.nlsr: Output model coefficients for nlsr object.}
\usage{
\method{coef}{nlsr}(object, ...)
}
\arguments{
\item{object}{An object of class 'nlsr'}

\item{...}{dot-args to provide exogenous data to the problem
Any data needed for the function. We do not know of any!}
}
\value{
returns the coefficients from the nlsr object. ?? as a named vector
}
\description{
A routine to extract and display the coefficients for a model 
estimated by \code{nlxb} or \code{nlfb} in the \code{nlsr} structured
\code{object}.
}
\section{Usage}{


   coef(object, ...)

\arguments{
\item{object}{An object of class 'nlsr'}
\item{\dots}{Any data needed for the function. We do not know of any!}
}
}

\section{Details}{


\code{coef.nlsr} extracts and displays the coefficients for a model 
    estimated by \code{nlxb} or \code{nlfb}.
}

\seealso{
Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
   
\keyword{nonlinear least squares}
}
\author{
John C Nash <nashjc@uottawa.ca>
}
