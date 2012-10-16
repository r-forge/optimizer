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
\alias{compl.par}
\alias{print.compl.par}
\alias{summary.compl.par}
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

phiMan(a, b, f, fprime)
GrAphiMan(a, b, f, fprime)
GrBphiMan(a, b, f, fprime)

phiKK(a, b, lambda)
GrAphiKK(a, b, lambda)
GrBphiKK(a, b, lambda)


phiLT(a, b, q)
GrAphiLT(a, b, q)
GrBphiLT(a, b, q)

compl.par(type=c("FB", "Min", "Man", "LT", "KK"), 
f, fprime, q, lambda)

\method{print}{compl.par}(x, \dots)

\method{summary}{compl.par}(object, \dots)



}
\arguments{
  \item{a}{first parameter.}
  \item{b}{second parameter.}
  \item{f, fprime}{a univariate function and its derivative.}
  \item{lambda}{a parameter in [0, 2[.}
  \item{q}{a parameter >1.}
  
  \item{type}{a character string for the complementarity
	function type: either \code{"FB"}, \code{"Min"}, \code{"Man"}, 
	\code{"LT"} or \code{"KK"}.}
  \item{x,object}{an object of class \code{"compl.par"}.}		
  \item{\dots}{further arguments to pass to \code{print}, 
		or \code{summary}.}			  
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

\code{compl.par} creates an object of class \code{"compl.par"} with attributes
\code{"type"} a character string and \code{"fun","grA","grB"} the corresponding 
functions for a given type. 
Optional arguments are also available, e.g. lambda for the KK complementarity 
function.
	
}
\value{
A numeric or an object of class \code{"compl.par"}.
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

complFB <- compl.par("FB") 
summary(complFB)

complKK <- compl.par("KK", lambda=1) 
summary(complKK)

complKK$fun(1, 1, complKK$lambda)
complFB$fun(1, 1)

}
\keyword{nonlinear}
\keyword{optimize}

