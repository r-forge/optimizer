\name{ufn}
\alias{ufn}
\encoding{UTF-8}
\title{Wrap user objective function for optimization tools}
\concept{minimization}
\concept{maximization}
\description{
        Provides a wrapper around user functions for nonlinear optimization
	to try to control for inadmissible arguments to user objective, gradient
	or hessian functions, as well as provide for scaling and maximization.
}
\usage{
ufn(par, fnuser, ps, fs=1.0, maximize=FALSE, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
 \item{ps}{A vector of scalings for the parameters. The REAL parameters are 
         the internal parameters multiplied by the scalings. The optimizers 
         use upar/ps internally, where upar are the user space parameters. 
         Note that there MUST be a vector ps. Its length and properties are 
         NOT checked inside ufn, ugr or uhess for efficiency.}
 \item{fs}{A user supplied scaling for the function. Output is the fnuser$fn/fs.}
 \item{maximize}{Set TRUE if the objective function is to be maximized. Default FALSE.}
 \item{...}{Other data needed to evaluate the user function.}
}
\details{
  None
}
\value{
  \code{ufn} returns a scalar numeric value, but this is set to the R constant
   .Machine$double.xmax if the inputs to the function are inadmissible and the
   computation of \code{fn} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}
\examples{
cat("Show how ufn traps an inadmissible set of parameters to a user function\n")

badlogf<-function(x, skale=10){
   sq<-seq(1:length(x))
   r<-(10-x)^2 + skale*log(x-sq)
   f<-as.double(crossprod(r))
} # note that this will fail when length(x)>x for some element of x

badlogg<-function(x, skale=10){# This is the gradient of badlogf
   sq<-seq(1:length(x))
   r<-(10-x)^2 + skale*log(x-sq)
   g<-2*r*(-2*(10-x)+skale/(x-sq))
} # note that this will fail when length(x)>x for some element of x

#badlogh<-function(x, skale=10){
#   sq<-seq(1:length(x))
#   r<-(10-x)^2 + skale*log(x-sq)
#   H<-r%*%t(r) # WRONG!
#   2*r*(-2*(10-x)+skale/(x-sq))
## NOT YET SET UP PROPERLY #  
#} # note that this will fail when length(x)>x for some element of x

##badlog<-list(fn=badlogf, gr=badlogg, hess=badlogh)
badlog<-list(fn=badlogf, gr=badlogg, hess=NULL)

x0<-rep(20, 4)
ps1<-rep(1,4)
cat("OK parameters:")
print(x0)
fval<-ufn(x0, badlog, ps1)
print("result:")
print(fval)

skale=1
x0<-rep(20, 4)
cat("skale=",skale,"  OK parameters:")
print(x0)
fval<-ufn(x0, badlog, ps1, skale=skale)
print("result:")
print(fval)

x0<-rep(2, 4)
cat("Bad parameters:")
print(x0)
fval<-ufn(x0, badlog, ps1)
print("result:")
print(fval)

skale=1
x0<-rep(2, 4)
cat("skale=",skale,"  Bad parameters:")
print(x0)
fval<-ufn(x0, badlog, ps1, skale=skale)
print("result:")
print(fval)

genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
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
	gg
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

mygenrose<-list(fn=genrose.f, gr=genrose.g, hess<-genrose.h)

x0<-c(1,2,3,4)
x1<-c(1,1,1,1)
psc<-c(1,2,3,4)
ps1<-c(1,1,1,1)

g0<-genrose.f(x0)
g0
g1<-genrose.f(x1)
g1

tryf01<-ufn(x0, mygenrose, ps1, fs=1, gs=10)
tryf01
tryf11<-ufn(x1, mygenrose, ps1, fs=1, gs=10)
tryf11

tryf02<-ufn(x0, mygenrose, ps1, fs=2, gs=10)
tryf02
tryf12<-ufn(x1, mygenrose, ps1, fs=2, gs=10)
tryf12

tryf0c<-ufn(x0, mygenrose, psc, fs=1, gs=10)
tryf0c
tryf1c<-ufn(x1, mygenrose, psc, fs=1, gs=10)
tryf1c

tryf0c2<-ufn(x0, mygenrose, psc, fs=2, gs=10)
tryf0c2
tryf1c2<-ufn(x1, mygenrose, psc, fs=2, gs=10)
tryf1c2


# scaled quadratic function

par<-rep(pi,4)
sqf.f<-function(par,sbase=2){
   npar<-length(par)
   pshift<-seq(1,npar)
   ss<-sbase^(2*pshift)
   ff<-sum(ss*(par-pshift)^2)
}

sqf.g<-function(par,sbase=2){
   npar<-length(par)
   pshift<-seq(1,npar)
   ss<-sbase^(2*pshift)
   gg<-2*ss*(par-pshift)
}

npar<-length(par)
sbase<-2
aou<-optim(par,sqf.f, control=list(trace=1), sbase=sbase)
aou
ao<-optim(par,sqf.f, control=list(trace=1, parscale=psc, fnscale=1), sbase=sbase)
ao
aob<-optim(par,sqf.f, gr=sqf.g, method="BFGS", control=list(trace=1,parscale=psc, fnscale=1), sbase=sbase)
aob
aobu<-optim(par,sqf.f, gr=sqf.g, method="BFGS", control=list(trace=1), sbase=sbase)
aobu


}
\keyword{nonlinear}
\keyword{optimize}

