\name{potential.reduction}
\alias{potential.ce}
\alias{gradpotential.ce}
\alias{psi.ce}
\alias{gradpsi.ce}
\title{Potential reduction algorithm utility functions}
\description{
		Functions for the potential reduction algorithm
	}
\usage{

potential.ce(u, n, zeta)

gradpotential.ce(u, n, zeta)	

psi.ce(z, dimx, dimlam, Hfinal, argfun, zeta)

gradpsi.ce(z, dimx, dimlam, Hfinal, jacHfinal, argfun, argjac, zeta)


}
\arguments{
  \item{u}{a numeric vector : \eqn{u=(u_1, u_2)} 
	where \eqn{u_1} is of size \code{n}.}		  
  \item{n}{a numeric for the size of \eqn{u_1}.}
  \item{zeta}{a positive parameter.}
  
  \item{z}{a numeric vector : \eqn{z=(x, lambda, w)} 
	where \code{dimx} is the size of components of \eqn{x}
	and \code{dimlam} is the size of components of \eqn{lambda} and \eqn{w}.}
  \item{dimx}{a numeric vector with the size of each components of \eqn{x}.}	
  \item{dimlam}{a numeric vector with the size of each components of \eqn{lambda}.
	We must have \code{length(dimx) == length(dimlam)}.}
  \item{Hfinal}{the root function.}
  \item{argfun}{a list of additionnals arguments for \code{Hfinal}.}
  \item{jacHfinal}{the Jacobian of the root function.}
  \item{argjac}{a list of additionnals arguments for \code{jacHfinal}.}
}
\details{

\code{potential.ce} is the potential function for the GNEP, and \code{gradpotential.ce} its gradient.
\code{psi.ce} is the application of the potential function for \code{Hfinal}, and \code{gradpsi.ce}
its gradient.
	
}
\value{
A numeric or a numeric vector.
}
\references{

 S. Bellavia, M. Macconi, B. Morini (2003),
 \emph{An affine scaling trust-region approach to bound-constrained nonlinear systems},
 Applied Numerical Mathematics 44, 257-280

 A. Dreves, F. Facchinei, C. Kanzow and S. Sagratella (2011), 
 \emph{On the solutions of the KKT conditions of generalized Nash equilibrium problems}, 
 SIAM Journal on Optimization 21(3), 1082-1108.
 
} 
\seealso{
  See also \code{\link{GNE.ceq}}.
}
\author{
 Christophe Dutang
}

\keyword{nonlinear}
\keyword{optimize}

