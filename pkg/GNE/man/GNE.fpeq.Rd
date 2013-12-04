\name{GNE.fpeq}
\alias{GNE.fpeq}
\title{Fixed point equation reformulation of the GNE problem.}
\description{
		Fixed point equation reformulation via the NI function of the GNE problem.
	}
\usage{
GNE.fpeq(init, dimx, obj, argobj, grobj, arggrobj, 
	heobj, argheobj, joint, argjoint, jacjoint, argjacjoint, 
	method = "default", problem = c("NIR", "VIR"), 
	merit = c("NI", "VI", "FP"), order.method=1, control.outer=list(), 
	control.inner=list(), silent=TRUE, param=list(), stepfunc, argstep, ...)
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
	  
  \item{method}{either \code{"pure"}, \code{"UR"}, \code{"vH"}, \code{"RRE"}, \code{"MPE"}, 
	\code{"SqRRE"} or \code{"SqMPE"} method, see details. \code{"default"} 
	corresponds to \code{"MPE"}.}

  \item{problem}{either \code{"NIR"}, \code{"VIP"}, see details.}
  \item{merit}{either \code{"NI"}, \code{"VI"}, \code{"FP"}, see details.}
  \item{order.method}{the order of the extrapolation method.}
	
  \item{control.outer}{a list with control parameters for the fixed point algorithm.}
  \item{control.inner}{a list with control parameters for the fixed point function.}  
  \item{silent}{a logical to show some traces.}
  \item{param}{a list of parameters for the computation of the fixed point function.}
  
  \item{stepfunc}{the step function, only needed when \code{method="UR"}.}
  \item{argstep}{additional arguments for the step function.}
    
  \item{\dots}{further arguments to be passed to the optimization routine. 
	NOT to the functions.}
  
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


The fixed point approach consists in solving equation \eqn{y(x)=x}.

\describe{
	\item{(a) Crude or pure fixed point method:}{
		It simply consists in iterations \eqn{x_{n+1} = y(x_n)}.}
	\item{(b) Polynomial methods:}{
			\describe{
				\item{- relaxation algorithm (linear extrapolation):}{
					The next iterate is computed as \deqn{x_{n+1} = (1-\alpha_n) x_n + \alpha_n y(x_n).}
					The step \eqn{\alpha_n} can be computed in different ways: constant, decreasing
					serie or a line search method. In the literature of game theory, the decreasing serie
					refers to the method of Ursayev and Rubinstein (\code{method="UR"}) while the line search
					method refers to the method of von Heusinger (\code{method="vH"}). Note that the constant
					step can be done using the UR method.}
				\item{- RRE and MPE method:}{
					Reduced Rank Extrapolation and Minimal Polynomial Extrapolation 
					methods are polynomial extrapolation methods, where the monomials are functional 
					``powers'' of the y function, i.e. function composition of y. Of order 1, RRE and MPE
					consists of \deqn{x_{n+1} = x_n + t_n (y(x_n) - x_n),} 
					where \eqn{t_n} equals to
					\eqn{<v_n, r_n> / <v_n, v_n>} for RRE1 and \eqn{<r_n, r_n> / <v_n, r_n>} for MPE1, where
					\eqn{r_n =y(x_n) - x_n } and \eqn{v_n = y(y(x_n)) - 2y(x_n) + x_n}. 
					To use RRE/MPE methods, set \code{method = "RRE"} or \code{method = "MPE"}.}	
				\item{- squaring method:}{
					It consists in using an extrapolation method (such as RRE and MPE)
					after two iteration of the linear extrapolation, i.e. 
					\deqn{x_{n+1} = x_n -2 t_n r_n + t_n^2 v_n.} The squared version of RRE/MPE methods are
					available via setting \code{method = "SqRRE"} or \code{method = "SqMPE"}.}	
			}
		}
	\item{(c) Epsilon algorithms:}{Not implemented.}
}
For details on fixed point methods, see Varadhan & Roland (2004).

The \code{control.outer} argument is a list that can supply any of the following components:
\describe{
	\item{\code{merit="FP"} and \code{method="pure"}}{see \code{\link{fpiter}}.
		the default parameters are \code{list(tol=1e-6, maxiter=100, trace=TRUE)}.
		}
	\item{\code{merit="FP"} and \code{method!="pure"}}{see \code{\link{squarem}}.
		the default parameters are \code{list(tol=1e-6, maxiter=100, trace=TRUE)}.
		}
	\item{\code{merit!="FP"}}{parameters are 
		\describe{
			\item{\code{tol}}{The absolute convergence tolerance. Default to 1e-6.}
			\item{\code{maxit}}{The maximum number of iterations. Default to 100.}
			\item{\code{echo}}{A logical or an integer (0, 1, 2, 3) to print traces. 
				Default to \code{FALSE}, i.e. 0.}
			\item{\code{sigma, beta}}{parameters for von Heusinger algorithm. 
				Default to 9/10 and 1/2 respectively.}
			}
		}

}


		
}
\value{
A list with components:
	\describe{
		\item{\code{par}}{The best set of parameters found.}
		\item{\code{value}}{The value of the merit function.}
		\item{\code{outer.counts}}{A two-element integer vector giving the number of 
			calls to fixed-point and merit functions respectively.}					
		\item{\code{outer.iter}}{The outer iteration number.}
		\item{\code{code}}{
			         The values returned are
         \describe{
			\item{\code{1}}{Function criterion is near zero.
			Convergence of function values has been achieved.}
			\item{\code{4}}{Iteration limit \code{maxit} exceeded.}
			\item{\code{100}}{an error in the execution.}
			}
		}
		\item{\code{inner.iter}}{The iteration number when 
			computing the fixed-point function.}	
		\item{\code{inner.counts}}{A two-element integer 
			vector giving the number of calls to the gap function and its gradient 
			when computing the fixed-point function.}			
		\item{\code{message}}{a string describing the termination code}	
	}
}
\references{
 A. von Heusinger (2009),
 \emph{Numerical Methods for the Solution of the Generalized Nash Equilibrium Problem},
 Ph. D. Thesis.
 
 A. von Heusinger & J. Kanzow (2009),
 \emph{Optimization reformulations of the generalized Nash equilibrium problem using Nikaido-Isoda-type functions},
 Comput Optim Appl .
 
 S. Uryasev & R.Y. Rubinstein (1994), 
 \emph{On relaxation algorithms in computation of noncooperative equilibria}, 
 IEEE Transactions on Automatic Control.
 
 R. Varadhan & C. Roland (2004),
 \emph{Squared Extrapolation Methods (SQUAREM): A New Class of Simple and Efficient Numerical 
 Schemes for Accelerating the Convergence of the EM Algorithm},
 Johns Hopkins University, Dept. of Biostatistics Working Papers.
 
}

\seealso{

See \code{\link{GNE.ceq}}, \code{\link{GNE.minpb}} and \code{\link{GNE.nseq}} 
for other approaches.


}
\author{
 Christophe Dutang
}
\keyword{nonlinear}
\keyword{optimize}

