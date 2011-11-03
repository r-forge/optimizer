\name{uhess}
\alias{uhess}
\encoding{UTF-8}
\title{Wrapper for user Hessian function for optimization tools}
\concept{minimization}
\concept{maximization}
\description{
Provides a wrapper around user analytic Hessian function for nonlinear optimization
	to try to control for inadmissible arguments to user that function,
	as well as provide for maximization.
}
\usage{
uhess(par, fnuser, ps, fs=1.0, maximize=FALSE, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
 \item{ps}{A vector of scalings for the parameters. The REAL parameters are the internal 
         parameters multiplied by the scalings. The optimizers use upar/ps internally, 
         where upar are the user space parameters. Note that there MUST be a vector ps.
         Its length and properties are NOT checked inside ufn, ugr or uhess for efficiency.}
 \item{fs}{A user supplied scaling for the function. Output is the fnuser$fn/fs.}
 \item{maximize}{Set TRUE if the objective function is to be maximized. Default FALSE.}
 \item{...}{Other data needed to evaluate the user function.}
}
\details{
  None
}
\value{
  \code{uhess} returns a matrix numeric value, but if the inputs to the function
    are inadmissible, all elements are set to the R 
   constant .Machine$double.xmax and the
   computation of \code{fn} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}\examples{
cat("uhess: try an inadmissible set of parameters to a user function\n")

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
   H<-r \%*\% t(r) #  note that we must escape percent signs
} # note that this will fail when length(x)>x for some element of x

badlog<-list(fn=badlogf, gr=badlogg, hess=badlogh)
 
ps1<-rep(1,4)

x0<-rep(20, 4)
cat("OK parameters:")
print(x0)
hval<-uhess(x0, badlog, ps1)
print("result:")
print(hval)

skale=1
x0<-rep(20, 4)
cat("skale=",skale,"  OK parameters:")
print(x0)
hval<-uhess(x0, badlog, ps1, skale=skale)
print("result:")
print(hval)

x0<-rep(2, 4)
cat("Bad parameters:")
print(x0)
hval<-uhess(x0, badlog, ps1)
print("result:")
print(hval)

skale=1
x0<-rep(2, 4)
cat("skale=",skale,"  Bad parameters:")
print(x0)
hval<-uhess(x0, badlog, ps1, skale=skale)
print("result:")
print(hval)

}
\keyword{nonlinear}
\keyword{optimize}
