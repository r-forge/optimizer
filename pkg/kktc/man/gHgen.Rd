% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{gHgen}
\alias{gHgen}
\encoding{UTF-8}
\title{Generate gradient and Hessian for a function at given parameters.}
\description{
  The optimx general-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function creates
  a data frame of results. \code{gHgens} extracts elements from this 
  data frame.
}
\usage{
  gHgen(par, fn, gr=NULL, hess=NULL, ktrace=FALSE, ...) 
}
\arguments{
 \item{par}{Set of parameters, assumed to be at a minimum of the function \code{fn}.}
 \item{fn}{Name of the objective function.}
 \item{gr}{(Optional) function to compute the gradient of the objective function. If present,
     we use the Jacobian of the gradient as the Hessian and avoid one layer of numerical
     approximation to the Hessian.} 
 \item{hess}{(Optional) function to compute the Hessian of the objective function. This
     is rarely available, but is included for completeness.}
 \item{ktrace}{A logical flag which, if TRUE, monitors the progress of gHgen. Default FALSE.}
 \item{\dots}{Extra data needed to compute the function, gradient and Hessian.}
}

\details{
  None
}
\value{
  \code{gHgen} a list of two items, the gradient vector and the hessian matrix.
}
\examples{
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

parx<-rep(1,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgen(parx,genrose.f)
cat("results of gHgen for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgen(parx,genrose.f,genrose.g)
cat("results of gHgen for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")
parx<-rep(0.9,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgen(parx,genrose.f,ktrace=TRUE, gs=9.4)
cat("results of gHgen with gs=",9.4," for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgen(parx,genrose.f,genrose.g, ktrace=TRUE)
cat("results of gHgen for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
gst<-5
cat("\n\nTest with full calling sequence and gs=",gst,"\n")
gengall<-gHgen(parx,genrose.f,genrose.g,genrose.h, ktrace=TRUE,gs=gst)
print(gengall)
}
\keyword{nonlinear}
\keyword{optimize}
