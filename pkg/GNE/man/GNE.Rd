\name{GNE}
\alias{GNE}
\alias{GNE}
\title{GNE package}
\description{
		Generalized Nash Equilibrium computational methods.
	}
\details{

Computing generalized Nash Equilibrium can be done in three different approaches.
\describe{
	\item{(i) fixed point approach}{It consists in solving equation \eqn{y(x)=x}. }
	\item{(ii) gap function minimization}{It consists in minimizing a gap function \eqn{min V(x)}. }
	\item{(iii) extended KKT system}{It consists in solving the extended Karush-Kuhn-Tucker 
		(KKT) system \eqn{\Phi(z)=0}. }
}
The three different approaches are implemented in three distinct function \code{\link{fixedpoint}}, 
\code{\link{NewtonKKT}} and \code{\link{minGap}}. 

Additionnal utitilty functions are also available:
\code{\link{rejection}}, \code{\link{projector}}, \code{\link{stepfunc}},
 \code{\link{complfunc}} and \code{\link{templatePhi}}.


TODO:
\describe{
\item{-}{write a pdf vignette.}
\item{-}{write the headers for all files, in particuler C files.}
}

}
\references{

 F. Facchinei, A. Fischer & V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.

 A. von Heusinger (2009),
 \emph{Numerical Methods for the Solution of the Generalized Nash Equilibrium Problem},
 Ph. D. Thesis.
 
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .


}
\author{
 Christophe Dutang
}
\keyword{math}


