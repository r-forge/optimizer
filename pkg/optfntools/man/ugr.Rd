\name{ugr}
\alias{ugr}
\encoding{UTF-8}
\title{Wrapper for user gradient function for optimization tools}
\concept{minimization}
\concept{maximization}
\description{
        Provides a wrapper around user gradient function for nonlinear optimization
	to try to control for inadmissible arguments to user objective, gradient
	or hessian functions, as well as provide for maximization.
}
\usage{
ugr(par, fnuser, ps, fs=1.0, maximize=FALSE, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess, as well as a callpos that indicates the R stack 
         position (see details). fn generates the scalar numerical value of the
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
   The calling routine MUST have a variable \code{kgr} defined. This is incremented in each
   call to \code{ugr}. The stack position of the calling routine / environment is best set
   in the \code{fnuser} list by using  

     callpos = sys.nframe()

   as in the examples. 
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

badlog<-list(fn=badlogf, gr=badlogg, hess=badlogh, callpos=sys.nframe())

ps1<-rep(1,4)

x0<-rep(20, 4)
cat("OK parameters:")
print(x0)
gval<-ugr(x0, badlog, ps1)
print("result:")
print(gval)

skale=1
x0<-rep(20, 4)
cat("skale=",skale,"  OK parameters:")
print(x0)
gval<-ugr(x0, badlog, ps1, skale=skale)
print("result:")
print(gval)

x0<-rep(2, 4)
cat("Bad parameters:")
print(x0)
gval<-ugr(x0, badlog, ps1)
print("result:")
print(gval)

skale=1
x0<-rep(2, 4)
cat("skale=",skale,"  Bad parameters:")
print(x0)
gval<-ugr(x0, badlog, ps1, skale=skale)
print("result:")
print(gval)

cat("=============================================\n")
require("numDeriv")
fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}
kgr<-0
x0<-c(-1.2,1)
ga<-grr(x0)
cat("Rosenbrock at ")
print(x0)
cat("Analytic gradient:")
print(ga)
opxfn<-list(fn=fr, gr=NULL, hess=NULL, callpos=sys.nframe())
deriv.approx=NULL # ensure not set
cat("Default numerical approx:")
gn<-ugr(x0,fnuser=opxfn)
print(gn)

deriv.approx=grfwd 
cat("grfwd numerical approx:")
gn<-ugr(x0,fnuser=opxfn)
print(gn)

deriv.approx=grback 
cat("grback numerical approx:")
gn<-ugr(x0,fnuser=opxfn)
print(gn)

deriv.approx=grcentral
cat("grcentral numerical approx:")
gn<-ugr(x0,fnuser=opxfn)
print(gn)

deriv.approx=grnd # numDeriv
cat("numDeriv numerical approx:")
gn<-ugr(x0,fnuser=opxfn)
print(gn)

cat("called ugr ",kgr," times\n")

}

\keyword{nonlinear}
\keyword{optimize}
