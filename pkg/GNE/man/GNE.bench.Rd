\name{bench.GNE}
\alias{bench.GNE}
\alias{bench.GNE.nseq}
\alias{bench.GNE.ceq}
\alias{bench.GNE.fpeq}
\alias{bench.GNE.minpb}
\title{Benchmark function}
\description{
		Benchmark function to compare GNE computational methods.
	}
\usage{
bench.GNE.nseq(xinit, ..., echo=FALSE, control=list())
bench.GNE.ceq(xinit, ..., echo=FALSE, control=list())
bench.GNE.fpeq(xinit, ..., echo=FALSE, control.outer=list(), 
	control.inner=list())
bench.GNE.minpb(xinit, ..., echo=FALSE, control.outer=list(), 
	control.inner=list())
}
\arguments{
  \item{xinit}{a numeric vector for the initial point.}
  \item{\dots}{further arguments to be passed to \code{GNE.nseq}, 
    \code{GNE.ceq}, \code{GNE.fpeq} or \code{GNE.minpb}. 
	NOT to the functions \code{func1} and \code{func2}.}
  \item{echo}{a logical to get some traces of the benchmark computation.}
  \item{control, control.outer, control.inner}{a list with control 
  parameters to be passed to \code{GNE.xxx} function.}
}	
\details{
Computing generalized Nash Equilibrium can be done in three different approaches.
\describe{
	\item{(i) extended KKT system}{It consists in solving the non smooth extended Karush-Kuhn-Tucker 
		(KKT) system \eqn{\Phi(z)=0}. \code{func1} is \eqn{Phi} and \code{func2} is \eqn{Jac Phi}.}
	\item{(ii) fixed point approach}{It consists in solving equation \eqn{y(x)=x}. \code{func1} is \eqn{y} and 
		\code{func2} is ?}
	\item{(iii) gap function minimization}{It consists in minimizing a gap function \eqn{min V(x)}. 
		\code{func1} is \eqn{V} and \code{func2} is \eqn{Grad V}.}
}


}
\value{
For \code{GNE.bench.ceq} and \code{GNE.bench.nseq}, a \code{data.frame} 
is returned with columns:
	\describe{
		\item{\code{method}}{the name of the method.}	
		\item{\code{fctcall}}{the number of calls of the function.}	
		\item{\code{jaccall}}{the number of calls of the Jacobian.}	
		\item{\code{comptime}}{the computation time.}	
		\item{\code{normFx}}{the norm of the merit function at the final iterate.}	
		\item{\code{code}}{the exit code.}	
		\item{\code{localmethods}}{the name of the local method.}	
		\item{\code{globalmethods}}{the name of the globalization method.}	
		\item{\code{x}}{the final iterate.}	
	}
For \code{GNE.bench.minpb}, a \code{data.frame} 
is returned with columns:
	\describe{
		\item{\code{method}}{the name of the method.}	
		\item{\code{minfncall.outer}}{the number of calls of the merit function.}	
		\item{\code{grminfncall.outer}}{the number of calls of the 
			gradient of the merit function.}	
		\item{\code{gapfncall.inner}}{the number of calls of the gap function.}				
		\item{\code{grgapfncall.outer}}{the number of calls of the 
			gradient of the gap function.}	
		\item{\code{comptime}}{the computation time.}	
		\item{\code{normFx}}{the norm of the merit function at the final iterate.}	
		\item{\code{code}}{the exit code.}	
		\item{\code{x}}{the final iterate.}	
	}
For \code{GNE.bench.fpeq}, a \code{data.frame} 
is returned with columns:
	\describe{
		\item{\code{method}}{the name of the method.}	
		\item{\code{fpfncall.outer}}{the number of calls of the fixed-point function.}	
		\item{\code{merfncall.outer}}{the number of calls of the merit function.}	
		\item{\code{gapfncall.inner}}{the number of calls of the gap function.}				
		\item{\code{grgapfncall.outer}}{the number of calls of the 
			gradient of the gap function.}	
		\item{\code{comptime}}{the computation time.}	
		\item{\code{normFx}}{the norm of the merit function at the final iterate.}	
		\item{\code{code}}{the exit code.}	
		\item{\code{x}}{the final iterate.}	
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
\seealso{

See \code{\link{GNE.fpeq}}, \code{\link{GNE.minpb}}, \code{\link{GNE.ceq}} 
and \code{\link{GNE.nseq}} for other approaches.

}
\author{
 Christophe Dutang
}
\keyword{math}


