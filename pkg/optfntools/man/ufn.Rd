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
ufn(par, fnuser)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{fnuser}{A user-supplied function object that has three sub-functions
         fn, gr, and hess. fn generates the scalar numerical value of the
         objective function, gr its vector valued gradient (or is NULL) and
         hess a numerical matrix for the Hessian (or is NULL).}
}
\details{
   The usual dot arguments (...) are subsumed in fnuser$dots to save
   complexity in the function call. Note that we need to unlist()
   these is the call to the actual user function.
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
## source("/home/work/R-optimtest/develmake/optfntools/R/ufn.R")
## source("/home/work/R-optimtest/develmake/optfntools/R/optstart.R")

badlogf<-function(x, skale=10){
   cat("in badlogf, skale=",skale,"\n")
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


x0<-rep(20, 4)
setup<-optstart(length(x0)) # Need to sort out how to set this up appropriately
cat("setup =")
print(str(setup))
print(ls(setup))

opxfn<-list(fn=badlogf, gr=badlogg, hess=NULL, OPCON=setup, dots=NULL)

ps1<-rep(1,4)
cat("skale= NULL, parameters:")
print(x0)
cat("Calling original function:")
fval0<-badlogf(x0)
cat(fval0,"\n")
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",setup$KFN,"\n")

cat("======================================\n")
tmp<-readline("set skale=1 AND reset opxfn")
skale<-1
opxfn<-list(fn=badlogf, gr=badlogg, hess=NULL, OPCON=setup, dots=list(skale=skale))
cat("prior to call, opxfn$dots:")
print(opxfn$dots)

cat("skale=",skale," parameters:")
print(x0)
fval0<-badlogf(x0, skale=skale)
cat(fval0,"\n")
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",setup$KFN,"\n")

cat("======================================\n")
tmp<-readline("Compare good and bad params")

x0<-rep(20, 4)
cat("skale=",skale,"  OK parameters:")
print(x0)
tfval0<-try(badlogf(x0))
print(tfval0)


opxfn<-list(fn=badlogf, gr=badlogg, hess=NULL, OPCON=setup, dots=NULL)
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",setup$KFN,"\n")

cat("======================================\n")
tmp<-readline("now bad params")
x0<-rep(2, 4)
cat("Bad parameters:")
opxfn<-list(fn=badlogf, gr=badlogg, hess=NULL, OPCON=setup, dots=NULL)
print(x0)
tfval0<-try(badlogf(x0))
print(tfval0)
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",setup$KFN,"\n")

cat("======================================\n")
rm(setup) # Try to remove the scratchpad

}
\keyword{nonlinear}
\keyword{optimize}

