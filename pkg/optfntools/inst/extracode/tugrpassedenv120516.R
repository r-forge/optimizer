rm(list=ls())
require(optfntools)
cat("Show how ugr works\n")
source("/home/work/R-optimtest/xdevel/optfntools/R/ugr.R")
source("/home/work/R-optimtest/xdevel/optfntools/R/ufn.R")
## source("/home/work/R-optimtest/develmake/optfntools/R/optstart.R")

cat("matrix function\n")

aa<-matrix(c(2,1,1,2),nrow=2)

myxp<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   f<-as.numeric((t(par) %*% A) %*% par)+(as.numeric(crossprod(par))-1)^2
}

myxpg<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   gg<-2.0*as.vector(A %*% par)+4.0*(as.numeric(crossprod(par))-1)*par
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
#??????? tested to here
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
#   H<-r%*%t(r) # WRONG!
#   2*r*(-2*(10-x)+skale/(x-sq))
## NOT YET SET UP PROPERLY #  
#} # note that this will fail when length(x)>x for some element of x


x0<-rep(20, 4)
npar<-4
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))

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
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))

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
tmp<-readline("set skale=1 AND reset opxfn")

# ??? fix up below here -- fn not gr

skale<-1
npar<-length(x0)
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(skale=skale)))
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


opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))
fval<-ufn(x0, opxfn)
print("result:")
print(fval)
cat("counter: kfn=",opxfn$KFN,"\n")

cat("======================================\n")
tmp<-readline("now bad params")
x0<-rep(2, 4)
npar<-length(x0)
cat("Bad parameters:")
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))
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

