\name{compl}
\alias{complfunc}
\alias{phiFB}
\alias{phiMin}
\alias{GrAphiFB}
\alias{GrBphiFB}
\alias{GrAphiMin}
\alias{GrBphiMin}
\title{Complementarity functions}
\description{
		Classic Complementarity functions
	}
\usage{
phiFB(a, b) 
GrAphiFB(a, b) 
GrBphiFB(a, b) 

phiMin(a, b) 
GrAphiMin(a, b) 
GrBphiMin(a, b)
}
\arguments{
  \item{a}{first parameter.}
  \item{b}{second parameter.}
}
\details{
\code{phiFB} corresponds to the Fischer-Burmeister complementarity function 
\eqn{\sqrt{a^2+b^2} - (a+b)}{sqrt(a^2+b^2) - (a+b)} and 
\code{phiMin} to the minimum Complementarity function \eqn{\min(a,b)}{min(a,b)}.

\code{GrAXXX} and \code{GrBXXX} implements the derivative of the complementarity 
function \code{XXX} with respect to \eqn{a} and \eqn{b} respectively.
	
}
\value{
A numeric.
}
\references{
	F. Facchinei and J.S. Pang,
	\emph{Finite-Dimensional Variational Inequalities and Complementarity Problems}, 
	Springer-Verlag (New York 2003).
} 
\seealso{
  See also \code{\link{NewtonKKT}}.
}
\author{
 Christophe Dutang
}
\keyword{nonlinear}
\keyword{optimize}

