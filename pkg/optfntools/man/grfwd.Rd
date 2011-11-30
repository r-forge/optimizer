\name{grfwd}
\alias{grfwd}

\title{Forward difference numerical gradient approximation.}

\description{
   \code{grfwd} computes the forward difference approximation to the gradient of 
          user function \code{userfn}. 
}


\usage{
   grfwd(par, userfn, fbase=NULL, eps=1.0e-7, ...)
}


\arguments{
     
    \item{par}{
        parameters to the user objective function userfn
        }
    \item{userfn}{
        User-supplied objective function
        }
    \item{fbase}{
        The value of the function at the parameters, else NULL. This is to save
        recomputing the function at this point.
        }
    \item{eps}{
        a small number used to generate the finite difference step
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}


\details{
\tabular{ll}{
Package: \tab grfwd\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
}


\value{
\code{grfwd} returns a single vector object \code{df} which approximates the 
gradient of userfn at the parameters par.

}


\author{

    John C. Nash
}

\examples{
# Want examples and comparisons with numDeriv and other forms. ??

require("numDeriv")
cat("testgrprogs.R\n")

# genrosa function code -- attempts to match the rosenbrock at gs=100 and x=c(-1.2,1)
genrosa.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
        # Note do not at 1.0 so min at 0
	fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

genrosa.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn1]
        # f = gs*z1*z1 + z2*z2
	gg[tn] <- 2 * (gs * z1)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 - 2 *z2 
	return(gg)
}

genrosa.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
#		z2<-1.0 - x[i-1]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}


maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

gmaxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	g<- -4*(x-ss)*as.numeric(crossprod(x-ss))
	return(g)
}


negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

gnegmaxfn<-function(x) {
	g<-(-1)*gmaxfn(x)
	return(g)
}



cat("Genrosa test at ")

#?? could vary eps too??
x0<-rep(pi,4)
print(x0)

f0<-genrosa.f(x0)
cat("f0=",f0,"\n")
ga<-genrosa.g(x0)
cat("analytic gradient:")
print(ga)

gnnd<-grad(genrosa.f, x0)
cat("numDeriv gradient:")
print(gnnd)
cat("Max abs deviation: ",max(abs(ga-gnnd)),"\n")



gnfwd<-grfwd(x0,genrosa.f,fbase=f0)
cat("grfwd gradient:")
print(gnfwd)
cat("Max abs deviation: ",max(abs(ga-gnfwd)),"\n")

gnback<-grback(x0,genrosa.f,fbase=f0)
cat("grback gradient:")
print(gnback)
cat("Max abs deviation: ",max(abs(ga-gnback)),"\n")

gnctl<-grcentral(x0,genrosa.f)
cat("grcentral gradient:")
print(gnctl)
cat("Max abs deviation: ",max(abs(ga-gnctl)),"\n")
cat("==========================================\n")
cat("Now test with maxfn")

x0<-rep(pi,4)
print(x0)

f0<-maxfn(x0)
cat("f0=",f0,"\n")
ga<-gmaxfn(x0)
cat("analytic gradient:")
print(ga)

gnnd<-grad(maxfn, x0)
cat("numDeriv gradient:")
print(gnnd)
cat("Max abs deviation: ",max(abs(ga-gnnd)),"\n")

gnfwd<-grfwd(x0,maxfn,fbase=f0)
cat("grfwd gradient:")
print(gnfwd)
cat("Max abs deviation: ",max(abs(ga-gnfwd)),"\n")

gnback<-grback(x0,maxfn,fbase=f0)
cat("grback gradient:")
print(gnback)
cat("Max abs deviation: ",max(abs(ga-gnback)),"\n")

gnctl<-grcentral(x0,maxfn)
cat("grcentral gradient:")
print(gnctl)
cat("Max abs deviation: ",max(abs(ga-gnctl)),"\n")
cat("==========================================\n")
cat("Now test with negmaxfn")

x0<-rep(pi,4)
print(x0)

f0<-negmaxfn(x0)
cat("f0=",f0,"\n")
ga<-gnegmaxfn(x0)
cat("analytic gradient:")
print(ga)

gnnd<-grad(negmaxfn, x0)
cat("numDeriv gradient:")
print(gnnd)
cat("Max abs deviation: ",max(abs(ga-gnnd)),"\n")



gnfwd<-grfwd(x0,negmaxfn,fbase=f0)
cat("grfwd gradient:")
print(gnfwd)
cat("Max abs deviation: ",max(abs(ga-gnfwd)),"\n")

gnback<-grback(x0,negmaxfn,fbase=f0)
cat("grback gradient:")
print(gnback)
cat("Max abs deviation: ",max(abs(ga-gnback)),"\n")

gnctl<-grcentral(x0,negmaxfn)
cat("grcentral gradient:")
print(gnctl)
cat("Max abs deviation: ",max(abs(ga-gnctl)),"\n")





}

\keyword{optimize}

