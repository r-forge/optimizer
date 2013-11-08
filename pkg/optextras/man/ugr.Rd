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
ugr(par, fnuser)
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
  \code{ugr} returns a vector numeric value, but all elements are set to the R 
   constant .Machine$double.xmax if the inputs to the function are inadmissible and the
   computation of \code{gr} fails. The returned value has an attribute 
   \code{inadmissible} which is returned TRUE in this case, but otherwise
   is FALSE.
}
\examples{
cat("Show how ugr works\n")
cat("matrix function\n")
aa<-matrix(c(2,1,1,2),nrow=2)

myxp<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   f<-as.numeric((t(par) \%*\% A) \%*\% par)+(as.numeric(crossprod(par))-1)^2
}

myxpg<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   gg<-2.0*as.vector(A \%*\% par)+4.0*(as.numeric(crossprod(par))-1)*par
}
npar<-2
opxfn<-list2env(list(fn=myxp, gr=myxpg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(A=aa)))

x0<-c(1,1)
g0<-myxpg(x0, A=aa)
print(g0)

cat("using numDeriv on myxp\n")
gn<-grad(myxp, x0, A=aa)
print(gn)

cat("using ugr with myxpg\n")
g0u<-ugr(x0,opxfn)
print(g0u)
rm(opxfn)
cat("using ugr with numderiv \"grnd\"\n")
opxfn<-list2env(list(fn=myxp, gr="grnd", hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(A=aa)))

g0un<-ugr(x0,opxfn)
print(g0un)


tmp<-readline("next")
rm(opxfn)
cat("=====================================\n\n")

badlogf<-function(x, skale=10){
#   cat("in badlogf, skale=",skale,"\n")
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
cat("Calling analytical badlogg:")
gval0<-badlogg(x0)
print(gval0)
gvalu<-ugr(x0, opxfn)
cat("result:")
print(gvalu)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grfwd"
gvalunf<-ugr(x0, opxfn)
cat("result from grfwd:")
print(gvalunf)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grcentral"
gvalunc<-ugr(x0, opxfn)
cat("result from grcentral:")
print(gvalunc)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grback"
gvalunb<-ugr(x0, opxfn)
cat("result from grback:")
print(gvalunb)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grnd"
gvalund<-ugr(x0, opxfn)
cat("result from grnd:")
print(gvalund)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

cat("======================================\n")
tmp<-readline("change parameter scaling")

x0<-rep(20, 4)
npar<-4
opxargs<-list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE,
    PARSCALE=rep(1,npar), FNSCALE=1, KFN=0, KGR=0, KHESS=0, dots=NULL)
opxfn<-list2env(opxargs)

ps1<-1/(1:4)
opxfn$PARSCALE<-ps1
cat("skale= NULL, parameters:")
print(x0)
cat("parscale:")
print(ps1)
cat("and in opxfn:")
print(opxfn$PARSCALE)
cat("Calling badlogg function:\n")
gval0<-badlogg(x0)
print(gval0)
cat("grad on badlogf:\n")
print(grad(badlogf,x0))

x0s<-x0/ps1
cat("x0s:")
print(x0s)

gvalu<-ugr(x0s, opxfn)
cat("result of ugr:\n")
print(gvalu)
cat("rescaled:")
print(gvalu/ps1)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")


opxfn$gr<-"grfwd"
gvalunf<-ugr(x0s, opxfn)
cat("result from grfwd:")
print(gvalunf)
cat("rescaled:")
print(gvalunf/ps1)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grcentral"
gvalunc<-ugr(x0s, opxfn)
cat("result from grcentral:")
print(gvalunc)
cat("rescaled:")
print(gvalunc/ps1)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grback"
gvalunb<-ugr(x0s, opxfn)
cat("result from grback:")
print(gvalunb)
cat("rescaled:")
print(gvalunb/ps1)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

opxfn$gr<-"grnd"
gvalund<-ugr(x0s, opxfn)
cat("result from grnd:")
print(gvalund)
cat("rescaled:")
print(gvalund/ps1)
cat("counter: kfn=",opxfn$KFN," kgr=",opxfn$KGR,"\n")

cat("======================================\n")

}

\keyword{nonlinear}
\keyword{optimize}
