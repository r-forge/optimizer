\name{NIR}
\alias{gapNIR}
\alias{gradxgapNIR}
\alias{gradygapNIR}
\alias{fpNIR}
\title{Nikaido Isoda Reformulation}
\description{
		functions of the Nikaido Isoda Reformulation of the GNEP
	}
\usage{

gapNIR(x, y, dimx, obj, argobj, param=list(), echo=FALSE)
gradxgapNIR(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
gradygapNIR(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
fpNIR(x, dimx, obj, argobj, joint, argjoint,  
	grobj, arggrobj, jacjoint, argjacjoint, param=list(), 
	echo=FALSE, control=list(), yinit=NULL, optim.method="default")


}
\arguments{
  \item{x,y}{a numeric vector.}
  \item{dimx}{a vector of dimension for \code{x}.}

  \item{obj}{objective function (to be minimized), see details.}
  \item{argobj}{a list of additional arguments.}
  \item{grobj}{gradient of the objective function, see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}

  \item{joint}{joint function, see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{jacjoint}{gradient of the joint function, see details.}
  \item{argjacjoint}{a list of additional arguments of the joint Jacobian.}

  \item{param}{ a list of parameters.}
  \item{control}{a list with control parameters for the fixed point algorithm.}
  \item{yinit}{initial point when computing the fixed-point function.}
  \item{optim.method}{optimization method when computing the fixed-point function.}

  \item{echo}{a logical to show some traces.}    
}
\details{
\code{gapNIR} computes the Nikaido Isoda function of the GNEP, while \code{gradxgapNIR} 
and \code{gradygapNIR} give its gradient with respect to \eqn{x} and \eqn{y}.
\code{fpNIR} computes the fixed-point function.
	
}
\value{
A vector for \code{funSSR} or a matrix for \code{jacSSR}.
}
\references{
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .
 
 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
} 
\seealso{
  See also \code{\link{GNE.fpeq}}.
}
\author{
 Christophe Dutang
}
\examples{

}
\keyword{math}
\keyword{optimize}

