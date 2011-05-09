\name{templatePhi}
\alias{templatePhi}
\alias{Phi}
\alias{phi}
\alias{JacPhi}
\alias{jacphi}
\alias{Jacphi}
\alias{jacPhi}
\title{the Phi function}
\description{
		a template for the Phi function of the extended KKT system
	}
\usage{
Phi(z, dimx, dimlam, grobj, arggrobj, constr, argconstr, grconstr, arggrconstr, compl)
JacPhi(z, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
	heconstr, argheconstr, gcompla, gcomplb)
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
  \item{compl}{the complementarity function with two arguments: \code{compl(a,b)}.}

  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}	
  \item{heconstr}{Hessian of the constraint function, see details.}
  \item{argheconstr}{a list of additional arguments of the constraint Hessian.}
  \item{gcompla}{derivative of the complementarity function w.r.t. the first argument.}
  \item{gcomplb}{derivative of the complementarity function w.r.t. the second argument.}    
}
\details{
Compute the Phi function or the Jacobian of the Phi function defined in F. Facchinei et al.(2009).

\describe{
    \item{\bold{Arguments of the Phi function}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{grobj}}{It must have 4 arguments: vector \code{z}, player number, derivative index and additional arguments \code{arggrobj}.}
		\item{\code{constr}}{It must have 4 arguments: vector \code{z}, player number, derivative index and additional arguments \code{argconstr}.}
		\item{\code{grconstr}}{It must have 4 arguments: vector \code{z}, player number, derivative index and additional arguments \code{arggrconstr}.}
		\item{\code{compl}}{It must have two arguments. A typical example is the minimum function.}
	}
	}
	\item{\bold{Arguments of the Jacobian of Phi}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{heobj}}{It must have 5 arguments: vector \code{z}, player number, two derivative indexes and additional arguments \code{argheobj}.}
		\item{\code{heconstr}}{It must have 5 arguments: vector \code{z}, player number, two derivative indexes and additional arguments \code{argheconstr}.}
		\item{\code{gcompla}}{It must have two arguments.}
	}
	}
}
See the example	below.
	
}
\value{
A vector for \code{Phi} or a matrix for \code{JacPhi}.
}
\references{
 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
} 
\seealso{
  See also \code{\link{NewtonKKT}}.
}
\author{
 Christophe Dutang
}
\examples{

# (1) objective linked functions
#

dimx <- c(2, 2, 3)

#Gr_x_j O_i(x)
grfullob <- function(x, i, j, par)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- 3*(x - 1:7)^2
	}
	if(i == 2)
	{
		grad <- 1:7*(x - 1:7)^(0:6)
	}
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2 - 5	
		grad <- c(1, 0, 1, 0, 4*x[5]*s, 4*x[6]*s, 4*x[7]*s)
			
	}
	grad[j]	+ par$a	
}


#Gr_x_k Gr_x_j O_i(x)
hefullob <- function(x, i, j, k, par)
{
	x <- x[1:7]
	if(i == 1)
	{
		he <- diag( 6*(x - 1:7) )
	}
	if(i == 2)
	{
		he <- diag( c(0, 2, 6, 12, 20, 30, 42)*(x - 1:7)^c(0, 0:5) )
	}
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2	
		
		he <- rbind(rep(0, 7), rep(0, 7), rep(0, 7), rep(0, 7),
			c(0, 0, 0, 0, 4*s+8*x[5]^2, 8*x[5]*x[6], 8*x[5]*x[7]),
			c(0, 0, 0, 0, 8*x[5]*x[6], 4*s+8*x[6]^2, 8*x[6]*x[7]),
			c(0, 0, 0, 0,  8*x[5]*x[7], 8*x[6]*x[7], 4*s+8*x[7]^2) )
	}
	he[j,k]	+ par$a	
}



# (2) constraint linked functions
#

dimlam <- c(1, 2, 2)

#constraint function g_i(x)
g <- function(x, i, par)
{
	x <- x[1:7]
	#cat(x[1:5], "|", i, "\n")
	if(i == 1)
		res <- sum( x^(1:7) ) -7
	if(i == 2)
		res <- c(sum(x) + prod(x) - 14, 20 - sum(x))
	if(i == 3)
		res <- c(sum(x^2) + 1, 100 - sum(x))
	#cat("res", res + par$a, "\n")	
	res	+ par$a				
}


#Gr_x_j g_i(x)
grfullg <- function(x, i, j, par)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- (1:7) * x ^ (0:6)
	}
	if(i == 2)
	{
		grad <- 1 + sapply(1:7, function(i) prod(x[-i]))
		grad <- cbind(grad, -1)
	}
	if(i == 3)
	{
		grad <- cbind(2*x, -1)
	}


	if(i == 1)
		res <- grad[j]	+ par$a				
	if(i != 1)
		res <- grad[j,]	+ par$a		
	as.numeric(res)
}



#Gr_x_k Gr_x_j g_i(x)
hefullg <- function(x, i, j, k, par)
{
	x <- x[1:7]
	if(i == 1)
	{
		he1 <- diag( c(0, 2, 6, 12, 20, 30, 42) * x ^ c(0, 0, 1:5) )
	}
	if(i == 2)
	{
		he1 <- matrix(0, 7, 7)
		he1[1, -1] <- sapply(2:7, function(i) prod(x[-c(1, i)]))
		he1[2, -2] <- sapply(c(1, 3:7), function(i) prod(x[-c(2, i)]))
		he1[3, -3] <- sapply(c(1:2, 4:7), function(i) prod(x[-c(3, i)]))
		he1[4, -4] <- sapply(c(1:3, 5:7), function(i) prod(x[-c(4, i)]))
		he1[5, -5] <- sapply(c(1:4, 6:7), function(i) prod(x[-c(5, i)]))
		he1[6, -6] <- sapply(c(1:5, 7:7), function(i) prod(x[-c(6, i)]))
		he1[7, -7] <- sapply(1:6, function(i) prod(x[-c(7, i)]))
						
						
		he2 <- matrix(0, 7, 7)
		
	}
	if(i == 3)
	{
		he1 <- diag(rep(2, 7))
		he2 <- matrix(0, 7, 7)
	}
	if(i != 1)
		return( c(he1[j, k], he2[j, k])	+ par$a	)
	else				
		return( he1[j, k]	+ par$a	)
}


FB <- function(a,b) sqrt(a^2+b^2)-(a+b)
GrAphiFB <- function(a, b) a / sqrt(a^2+b^2) - 1
GrBphiFB <- function(a, b) b / sqrt(a^2+b^2) - 1


# (3) compute Phi
#

z <- rexp(sum(dimx) + sum(dimlam))

a <- 0

Phi(z, dimx, dimlam, grfullob, list(a=a),  
	g, list(a=a), grfullg, list(a=a), FB)
	
# (4) compute Jac Phi
#
	
JacPhi(z, dimx, dimlam, hefullob, list(a=a),  
	g, list(a=a), grfullg, list(a=a),
	hefullg, list(a=a), GrAphiFB, GrBphiFB)

}
\keyword{math}
\keyword{optimize}

