% File src/library/stats/man/kktc.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{kktc}
\alias{kktc}
\encoding{UTF-8}
\title{Check Kuhn Karush Tucker conditions for a supposed function minimum}
\concept{minimization}
\concept{maximization}
\description{
  Provide a check on Kuhn-Karush-Tucker conditions based on quantities
  already computed. Some of these used only for reporting.
}
\usage{
     kktc(par, fval, ngmax, nHes, 
             control=list(kkttol=1e-3, kkt2tol=1e-6, ktrace=FALSE) )
}
\arguments{
 \item{par}{A vector of values for the parameters which are supposedly optimal.}
 \item{fval}{Function value at parameters \code{par}.}
 \item{ngmax}{Maximum absolute component of the gradient of the function at 
       \code{par}.}
 \item{nHes}{Matrix that is supposedly the Hessian at \code{par}.}
 \item{control}{A list of controls for the function, made up of the following items}
    \itemize{
       \item{kkttol}{A tolerance for the KKT gradient check.}
       \item{kkt2tol}{A tolerance used for the KKT Hessian check.}
       \item{ktrace}{A logical flag that, if TRUE, displays intermediate information
          on the progress of the function. Default is FALSE.}
    }
}
\details{
  At the moment, kktc does NOT handle box constraints or masks, thought the test built into
the \code{optimix} package does do this. The present code was developed by lifting it from \code{optimix}
so KKT checks could be added to other optimization routines.
}

\value{
  The output is a list consisting of 
  \item{ngmax}{The absolute value of the largest gradient component in magnitude.}
  \item{evratio}{The ratio of the smallest to largest Hessian eigenvalue. Note that this
     may be negative.}
  \item{kkt1}{A logical value that is TRUE if we consider the first (i.e., gradient) 
     KKT condition to be satisfied. WARNING: The decision is dependent on tolerances and
     scaling that may be inappropriate for some problems.}
  \item{kkt2}{A logical value that is TRUE if we consider the second (i.e., positive
      definite Hessian) KKT condition to be satisfied. WARNING: The decision is dependent 
      on tolerances and scaling that may be inappropriate for some problems.}
}
\seealso{

    \code{\link{optim}}

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
out<-kktc(parx,fval,max(abs(gval)),Ahess, control=list(ktrace=TRUE))
cat("results of kktc for genrose at ")
print(parx)
print(out)
cat("\n\n")
parx<-rep(0.9,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
out<-kktc(parx,fval,max(abs(gval)),Ahess, control=list(ktrace=TRUE))
cat("results of kktc for genrose at ")
print(parx)
print(out)
cat("\n\n")
cat("Example with optimization\n")
gst<-3
ans<-optim(parx,genrose.f,gs=gst)
print(ans)
paro<-ans$par
fval<-ans$value
gval<-grad(genrose.f,paro,gs=gst)
Ahess<-hessian(genrose.f,paro,gs=gst)
out<-kktc(parx,fval,max(abs(gval)),Ahess, control=list(ktrace=TRUE))
cat("results of kktc for genrose at ")
print(paro)
print(out)

}
\keyword{nonlinear}
\keyword{optimize}
