\name{eqsolve}
\alias{eqsolve}
\title{Solving non linear equations}
\description{
		Non linear Solving methods
	}
\usage{
eqsolve(xinit, f, jac, method=c("Newton", "Levenberg-Marquardt","Broyden"), 
	global=c("line search", "trust region", "none"), control=list())
}
\arguments{
  \item{xinit}{initial point.}
  \item{f}{the function for which we search roots.}
  \item{jac}{the Jacobian of the function \code{f}.}
  \item{method}{a character string specifying the method to use: either 
  \code{"Newton"}, \code{"Levenberg-Marquardt"}, or \code{"Broyden"}.}	
  \item{global}{a character string for the globalization method to be used:
  either \code{"line search"}, \code{"trust region"}, or \code{"none"}.}
  \item{control}{a list for the control parameters. See details.}
}
\details{

The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
	\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
	\item{\code{echo}}{A logical or an integer (0, 1, 2, 3, 4) to print traces. 
		Default to \code{FALSE}, i.e. 0.}
	\item{\code{echofile}}{A character string to store the traces in that file. Default to \code{NULL}.}	
	\item{\code{echograph}}{A character string to plot iter-by-iter information. 
		Either \code{"NULL"} (default), or \code{"line"} for line search plot 
		or \code{"trust"} for trust region plots.}	
	\item{\code{sigma}}{Reduction factor for the geometric linesearch. Default to 0.5.}
	\item{\code{btol}}{The backtracking tolerance. Default to 0.01.}	
	\item{\code{delta}}{The exponent parameter for the LM parameter, should in \eqn{[1,2]}. Default to 2.}
	\item{\code{initlnsrch}}{The initial integer for starting the line search. Default to 0.}		
	\item{\code{minstep}}{The minimal step. Default to 0.001.}		
}

}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{counts}}{A two-element integer vector giving the number of calls to \code{phi} and \code{jacphi} respectively.}			
		\item{\code{iter}}{The iteration number.}
		\item{\code{code}}{0 if convergence, 1 if \code{maxit} is reached, 10 
			if \code{tol} is not reached and 11 for both.}
	}
}

\seealso{
See \code{\link{nleqslv}} from the package of the same name.
}
\author{
 Christophe Dutang
}
\examples{

a <- 1

}
\keyword{nonlinear}
\keyword{optimize}

