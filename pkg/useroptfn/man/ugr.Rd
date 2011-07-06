\name{ugr}
\alias{ugr}
\encoding{UTF-8}
\title{Output tool for answers from optimx}
\concept{minimization}
\concept{maximization}
\description{
Provides a wrapper around user functions for nonlinear optimization
	to try to control for inadmissible arguments to user objective, gradient
	or hessian functions, as well as provide for maximization.
}
\usage{
ugr(par, fnuser, maximize=FALSE, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
 \item{maximize}{Set TRUE if the objective function is to be maximized. Default FALSE.}
 \item{...}{Other data needed to evaluate the user function.}
}
\details{
  None
}
\value{
  \code{ugr} returns a vector numeric value, but all elements are set to the R 
   constant .Machine$double.xmax if the inputs to the function are inadmissible and the
   computation of \code{gr} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}
\examples{
cat("Show how ugr traps an inadmissible set of parameters to a user function\n")
# currently no maximize??

badlogf<-function(x, skale=10){
   sq<-seq(1:length(x))
   r<-(10-x)^2 + skale*log(x-sq)
   f<-crossprod(r)
} # note that this will fail when length(x)>x for some element of x

badlogg<-function(x, skale=10){
   sq<-seq(1:length(x))
   r<-(10-x)^2 + skale*log(x-sq)
} # note that this will fail when length(x)>x for some element of x

badlogh<-function(x, skale=10){
   sq<-seq(1:length(x))
   r<-(10-x)^2 + skale*log(x-sq)
   H<-r%*%t(r)
} # note that this will fail when length(x)>x for some element of x

badlog<-list(fn=badlogf, gr=badlogg, hess=badlogh)

x0<-rep(20, 4)
cat("OK parameters:")
print(x0)
gval<-ugr(x0, badlog)
print("result:")
print(gval)

skale=1
x0<-rep(20, 4)
cat("skale=",skale,"  OK parameters:")
print(x0)
gval<-ugr(x0, badlog, skale=skale)
print("result:")
print(gval)

x0<-rep(2, 4)
cat("Bad parameters:")
print(x0)
gval<-ugr(x0, badlog)
print("result:")
print(gval)

skale=1
x0<-rep(2, 4)
cat("skale=",skale,"  Bad parameters:")
print(x0)
gval<-ugr(x0, badlog, skale=skale)
print("result:")
print(gval)

}

\keyword{nonlinear}
\keyword{optimize}
