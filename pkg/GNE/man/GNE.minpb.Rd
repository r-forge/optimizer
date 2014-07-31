\name{GNE.minpb}
\alias{GNE.minpb}
\title{Non smooth equation reformulation of the GNE problem.}
\description{
		Non smooth equation reformulation via the extended KKT system of the GNE problem.
	}
\usage{
GNE.minpb(init, dimx, obj, argobj, grobj, arggrobj, 
	heobj, argheobj, joint, argjoint, jacjoint, argjacjoint, 
	method="default", problem = c("NIR", "VIR"), control.outer=list(), 
	control.inner=list(), silent=TRUE, param=list(), 
	optim.type=c("free","constr"), ...)
}
\arguments{
  \item{init}{Initial values for the parameters to be optimized over: \eqn{z=(x, lambda, mu)}.}
  \item{dimx}{a vector of dimension for \eqn{x}.}
  
  \item{obj}{objective function (to be minimized), see details.}
  \item{argobj}{a list of additional arguments.}
  \item{grobj}{gradient of the objective function, see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}
  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}

  \item{joint}{joint function (\eqn{h(x)<=0}), see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{jacjoint}{Jacobian of the joint function, see details.}
  \item{argjacjoint}{a list of additional arguments of the Jacobian.}  

  \item{method}{either \code{"BB"}, \code{"CG"} or \code{"BFGS"}, see details.}
  \item{problem}{either \code{"NIR"}, \code{"VIP"}, see details.}
  \item{optim.type}{either \code{"free"}, \code{"constr"}, see details.}
  
  \item{control.outer}{a list with control parameters for the minimization algorithm.}
  \item{control.inner}{a list with control parameters for the minimization function.}  

  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions \code{phi} and \code{jacphi}.}

  \item{silent}{a logical to show some traces.}
  \item{param}{a list of parameters for the computation of the minimization function.}
  
}
\details{


Functions in argument must respect the following template:
\itemize{
  \item{\code{obj} must have arguments the current iterate \code{z}, the player number \code{i} 
  and optionnally additional arguments given in a list.}
  \item{\code{grobj} must have arguments the current iterate \code{z}, the player number \code{i}, 
  the derivative index \code{j} and optionnally additional arguments given in a list.}
  \item{\code{heobj} must have arguments the current iterate \code{z}, the player number \code{i}, 
  the derivative indexes \code{j}, \code{k} and optionnally additional arguments given in a list.}
  \item{\code{joint} must have arguments the current iterate \code{z}
  and optionnally additional arguments given in a list.}
  \item{\code{jacjoint} must have arguments the current iterate \code{z},
  the derivative index \code{j} and optionnally additional arguments given in a list.}
}



The gap function minimization consists in minimizing a gap function \eqn{min V(x)}. The function \code{minGap}
provides two optimization methods to solve this minimization problem.
\describe{
	\item{Barzilai-Borwein algorithm}{when \code{method = "BB"}, we use Barzilai-Borwein iterative scheme
	to find the minimum.}
	\item{Conjugate gradient algorithm}{when \code{method = "CG"}, we use the CG iterative
	scheme implemented in \code{R}, an Hessian-free method. }	
	\item{Broyden-Fletcher-Goldfarb-Shanno algorithm}{when \code{method = "BFGS"}, we use the BFGS iterative
	scheme implemented in \code{R}, a quasi-Newton method with line search. }	
}
In the game theory literature, there are two main gap functions: the regularized 
Nikaido-Isoda (NI) function and the regularized QVI gap function.
This correspond to \code{type="NI"} and \code{type="VI"}, respectively. 
See von Heusinger & Kanzow (2009) for details on the NI function and 
Kubota & Fukushima (2009) for the QVI regularized gap function.


The \code{control.outer} argument is a list that can supply any of the following components:
\describe{
	\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
	\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
	\item{\code{echo}}{A logical or an integer (0, 1, 2, 3) to print traces. 
		Default to \code{FALSE}, i.e. 0.}
	\item{\code{stepinit}}{Initial step size for the BB method (should be 
		small if gradient is ``big''). Default to 1.}
}


Note that the \code{Gap} function can return a numeric or a list with computation details. In the
latter case, the object return must be a list with the following components
\code{value}, \code{counts}, \code{iter}, see the example below.

}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{value}}{The value of the merit function.}
		\item{\code{outer.counts}}{A two-element integer vector giving the number of 
			calls to \code{Gap} and \code{gradGap} respectively.}				
		\item{\code{outer.iter}}{The outer iteration number.}
		\item{\code{code}}{
			         The values returned are
         \describe{
			\item{\code{1}}{Function criterion is near zero.
			Convergence of function values has been achieved.}
			\item{\code{2}}{x-values within tolerance. This means that the relative distance between two
			consecutive x-values is smaller than \code{xtol}.}
			\item{\code{3}}{No better point found.
			This means that the algorithm has stalled and cannot find an acceptable new point.
			This may or may not indicate acceptably small function values.}
			\item{\code{4}}{Iteration limit \code{maxit} exceeded.}
			\item{\code{5}}{Jacobian is too ill-conditioned.}
			\item{\code{6}}{Jacobian is singular.}
			\item{\code{100}}{an error in the execution.}
			}
		}
		\item{\code{inner.iter}}{The iteration number when 
			computing the minimization function.}	
		\item{\code{inner.counts}}{A two-element integer 
			vector giving the number of calls to the gap function and its gradient 
			when computing the minimization function.}			
		\item{\code{message}}{a string describing the termination code}	
	}
}
\references{
 A. von Heusinger (2009),
 \emph{Numerical Methods for the Solution of the Generalized Nash Equilibrium Problem},
 Ph. D. Thesis.
 
 A. von Heusinger and C. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .
 
 K. Kubota and M. Fukushima (2009), 
 \emph{Gap function approach to the generalized Nash Equilibrium problem},
 Journal of Optimization theory and applications.
}
\seealso{

See \code{\link{GNE.fpeq}}, \code{\link{GNE.ceq}} and \code{\link{GNE.nseq}} 
for other approaches.

}
\author{
 Christophe Dutang
}
\keyword{nonlinear}
\keyword{optimize}

