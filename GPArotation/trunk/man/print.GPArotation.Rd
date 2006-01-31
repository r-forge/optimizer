\name{print.GPArotation}
\alias{print.GPArotation}
\alias{summary.GPArotation}
\alias{print.summary.GPArotation}

\title{Print and Summary Methods for GPArotation}
\usage{
    \method{print}{GPArotation}(x, digits=3, Table=FALSE, ...)
    \method{summary}{GPArotation}(object, ...)
    \method{print}{summary.GPArotation}(x, digits=3, ...)
}
\arguments{
    \item{object}{a GPArotation object to summarize.}
    \item{x}{a summary.GPArotation to print.}
    \item{digits}{precision of printed numbers.}
    \item{...}{further arguments passed to other methods.}
}

\value{The object printed or a summary object.}
\description{
   Print an object or summary of an object returned by \code{GPForth} 
   or \code{GPFoblq}.
}
\seealso{
   \code{\link{GPForth}}, 
   \code{\link[base]{summary}}
}

\concept{rotation}
\keyword{internal}
