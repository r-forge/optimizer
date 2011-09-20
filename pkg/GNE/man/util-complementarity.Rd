\name{compl}
\alias{complementarity}
\alias{phiFB}
\alias{phiMin}
\alias{phiMan}
\alias{phiKK}
\alias{phiLT}
\alias{GrAphiFB}
\alias{GrBphiFB}
\alias{GrAphiMin}
\alias{GrBphiMin}
\alias{GrAphiMan}
\alias{GrBphiMan}
\alias{GrAphiKK}
\alias{GrBphiKK}
\alias{GrAphiLT}
\alias{GrBphiLT}
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

phiMan(a, b, f)
GrAphiMan(a, b, fprime)
GrBphiMan(a, b, fprime)

phiKK(a, b, lambda)
GrAphiKK(a, b, lambda)
GrBphiKK(a, b, lambda)


phiLT(a, b, q)
GrAphiLT(a, b, q)
GrBphiLT(a, b, q)

}
\arguments{
  \item{a}{first parameter.}
  \item{b}{second parameter.}
  \item{f, fprime}{a univariate function and its derivative.}
  \item{lambda}{a parameter in [0, 2[.}
  \item{q}{a parameter >1.}
}
\details{

We implement 5 complementarity functions From Facchinei & Pang (2003). 
\describe{
	\item{(i) \code{phiFB}}{ the Fischer-Burmeister complementarity function 
		\eqn{\sqrt{a^2+b^2} - (a+b)}{sqrt(a^2+b^2) - (a+b)}.}
	\item{(ii) \code{phiMin}}{ the minimum complementarity function \eqn{\min(a,b)}{min(a,b)}.}
	\item{(iii) \code{phiMan}}{ the Mangasarian's family of complementarity function \eqn{f(|a-b|) - f(a) - f(b)}, 
	typically \eqn{f(t)=t} or \eqn{f(t)=t^3}.}
	\item{(iv) \code{phiKK}}{ the Kanzow-Kleinmichel complementarity function 
		\eqn{(\sqrt( (a-b)^2 + 2*\lambda*a*b ) - (a+b) ) / (2-\lambda)}{(sqrt( (a-b)^2 + 2*lambda*a*b ) - (a+b) ) / (2-lambda)}.}
	\item{(v) \code{phiLT}}{ the Luo-Tseng complementarity function \eqn{(a^q + b^q)^(1/q) - (a+b)}.}
}				


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
  See also \code{\link{GNE.nseq}}.
}
\author{
 Christophe Dutang
}
\examples{

phiFB(1, 2)
phiLT(1, 2, 2)
phiKK(1, 2, 1)

-2*phiMin(1, 2)
phiMan(1, 2, function(t) t)


}
\keyword{nonlinear}
\keyword{optimize}

