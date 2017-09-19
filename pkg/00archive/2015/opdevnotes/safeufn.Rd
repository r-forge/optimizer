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
cat("matrix function\n")

aa<-matrix(c(2,1,1,2),nrow=2)

myxp<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   nn<-names(par)
   f<-as.numeric((t(par) \%*\% A) \%*\% par)+(as.numeric(crossprod(par))-1)^2
}

x0<-c(1,1)
f0<-myxp(x0, A=aa)
print(f0)

npar<-2
opxargs<-list(fn=myxp, gr=NULL, hess=NULL, MAXIMIZE=FALSE, 
       PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=list(A=aa))
opxfn<-list2env(opxargs)
axp1<-optim(x0, ufn, control=list(trace=2), fnuser=opxfn)
print(axp1)

cat("=====================================\n\n")
cat("Bad function -- fails when length(x)>x\n")

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
#   H<-r\%*\%t(r) # WRONG!
#   2*r*(-2*(10-x)+skale/(x-sq))
## NOT YET SET UP PROPERLY #  
#} # note that this will fail when length(x)>x for some element of x


x0<-rep(20, 4)
npar<-4
opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, 
    PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=NULL)
opxfn<-list2env(opxargs)

ps1<-rep(1,4)
cat("skale= NULL, parameters:")
print(x0)
cat("Calling original function:")
fval0<-badlogf(x0)
cat(fval0,"\n")
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
tmp<-readline("change parameter scaling")
x0<-rep(20, 4)
npar<-4
opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, 
    PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=NULL)
opxfn<-list2env(opxargs)

ps1<-1:4
opxfn$PARSCALE<-ps1
cat("skale= NULL, parameters:")
print(x0)
cat("Calling original function:")
fval0<-badlogf(x0*ps1)
cat(fval0,"\n")
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
tmp<-readline("set skale=1 AND reset opxfn")
skale<-1
npar<-length(x0)
opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE,
     PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=list(skale=skale))
opxfn<-list2env(opxargs)
cat("prior to call, opxfn$dots:")
print(opxfn$dots)

cat("skale=",skale," parameters:")
print(x0)
fval0<-badlogf(x0, skale=skale)
cat(fval0,"\n")
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
tmp<-readline("Compare good and bad params")

x0<-rep(20, 4)
npar<-length(x0)

cat("skale=",skale,"  OK parameters:")
print(x0)
tfval0<-try(badlogf(x0))
print(tfval0)

opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, 
   PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=NULL)
opxfn<-list2env(opxargs)
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
tmp<-readline("now bad params")
x0<-rep(2, 4)
npar<-length(x0)
cat("Bad parameters:")
opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, 
  PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=NULL)
opxfn<-list2env(opxargs)
print(x0)
tfval0<-try(badlogf(x0))
print(tfval0)
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
cat("Try to remove opxfn\n")
print(ls())
rm(opxfn) # Try to remove the scratchpad
print(ls())



}
\keyword{nonlinear}
\keyword{optimize}

