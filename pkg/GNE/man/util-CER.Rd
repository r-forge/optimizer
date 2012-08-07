\name{CER}
\alias{CER}
\alias{funCER}
\alias{jacCER}
\title{Constrained Equation Reformulation}
\description{
		functions of the Constrained Equation Reformulation of the GNEP
	}
\usage{

funCER(z, dimx, dimlam, 
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	echo=FALSE)

jacCER(z, dimx, dimlam,
	heobj, argheobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	hejoint, arghejoint,
	echo=FALSE)

}
\arguments{
  \item{z}{a numeric vector z containing x then lambda values.}
  \item{dimx}{dimension of x.}
  \item{dimlam}{dimension of lambda.}
  \item{grobj}{gradient of the objective function, see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}
  
  \item{constr}{constraint function, see details.}
  \item{argconstr}{a list of additional arguments of the constraint function.}
  \item{grconstr}{gradient of the constraint function, see details.}
  \item{arggrconstr}{a list of additional arguments of the constraint gradient.}
  
  \item{dimmu}{a vector of dimension for \eqn{mu}.}
  \item{joint}{joint function, see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{grjoint}{gradient of the joint function, see details.}
  \item{arggrjoint}{a list of additional arguments of the joint gradient.}
  
  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}	
  \item{heconstr}{Hessian of the constraint function, see details.}
  \item{argheconstr}{a list of additional arguments of the constraint Hessian.}
  \item{hejoint}{Hessian of the joint function, see details.}
  \item{arghejoint}{a list of additional arguments of the joint Hessian.}
  \item{echo}{a logical to show some traces.}    
  
}
\details{
Compute the H function or the Jacobian of the H function defined in Dreves et al.(2009).

\describe{
    \item{\bold{Arguments of the H function}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{grobj}}{
			The gradient \eqn{Grad Obj} of an objective function \eqn{Obj} (to be minimized) must have 3 arguments for \eqn{Grad Obj(z, playnum, ideriv)}: vector \code{z}, player number, derivative index
			, and optionnally additional arguments in \code{arggrobj}.
			}
		\item{\code{constr}}{
			The constraint function \eqn{g} must have 2 arguments: vector \code{z}, player number,
			such that \eqn{g(z, playnum) <= 0}. Optionnally, \eqn{g} may have additional arguments in \code{argconstr}.
			}
		\item{\code{grconstr}}{
			The gradient of the constraint function \eqn{g}  must have 3 arguments: vector \code{z}, player number, derivative index,
			and optionnally additional arguments in \code{arggrconstr}.
			}
	}
	}
	\item{\bold{Arguments of the Jacobian of H}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{heobj}}{It must have 4 arguments: vector \code{z}, player number, two derivative indexes.}
		\item{\code{heconstr}}{It must have 4 arguments: vector \code{z}, player number, two derivative indexes.}
	}
	Optionnally, \code{heobj} and \code{heconstr} can have additional arguments \code{argheobj} and \code{argheconstr}.
	}
}
See the example	below.
	
}
\value{
A vector for \code{funCER} or a matrix for \code{jacCER}.
}
\references{

 Dreves, A., Facchinei, F., Kanzow, C. and Sagratella, S. (2011), 
 \emph{On the solutions of the KKT conditions of generalized Nash equilibrium problems},
 SIAM Journal on Optimization.

 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
} 
\seealso{
  See also \code{\link{GNE.ceq}}.
}
\author{
 Christophe Dutang
}
\examples{



#-------------------------------------------------------------------------------
# (1) Example 5 of von Facchinei et al. (2007)
#-------------------------------------------------------------------------------

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j)
{
	if(i == 1)
		res <- c(2*(x[1]-1), 0)
	if(i == 2)
		res <- c(0, 2*(x[2]-1/2))
	res[j]	
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k)
	2 * (i == j && j == k)

dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	sum(x[1:2]) - 1
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	1
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0


x0 <- rep(0, sum(dimx))
z0 <- c(x0, 2, 2, max(10, 5-g(x0, 1) ), max(10, 5-g(x0, 2) ) )

#true value is (3/4, 1/4, 1/2, 1/2)
funCER(z0, dimx, dimlam, grobj=grobj, 
	constr=g, grconstr=grg)

jacCER(z0, dimx, dimlam, heobj=heobj, 
	constr=g, grconstr=grg, heconstr=heg)



#-------------------------------------------------------------------------------
# (2) Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------


#constants
myarg <- list(d= 20, lambda= 4, rho= 1)

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	res <- -arg$rho * x[i]
	if(i == j)
	res <- res + arg$d - arg$lambda - arg$rho*(x[1]+x[2])
	-res
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
	arg$rho * (i == j) + arg$rho * (j == k)	


dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	-x[i]
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	-1*(i == j)
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0

#true value is (16/3, 16/3, 0, 0) 

x0 <- rep(0, sum(dimx))
z0 <- c(x0, 2, 2, max(10, 5-g(x0, 1) ), max(10, 5-g(x0, 2) ) )


funCER(z0, dimx, dimlam, grobj=grobj, arggrobj=myarg, 
	constr=g, grconstr=grg)

jacCER(z0, dimx, dimlam, heobj=heobj, 
	argheobj=myarg, constr=g, grconstr=grg, heconstr=heg)





}
\keyword{math}
\keyword{optimize}

