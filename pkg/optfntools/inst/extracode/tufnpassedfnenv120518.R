rm(list=ls())
cat("Show how ufn traps an inadmissible set of parameters to a user function\n")
source("/home/work/R-optimtest/xdevel/optfntools/R/ufn.R")
## source("/home/work/R-optimtest/develmake/optfntools/R/optstart.R")

##?? check MAXIMIZE and FNSCALE



cat("matrix function\n")

aa<-matrix(c(2,1,1,2),nrow=2)

myxp<-function(par, A=NULL){
   if(is.null(A))stop("MUST have matrix A")
   nn<-names(par)
   f<-as.numeric((t(par) %*% A) %*% par)+(as.numeric(crossprod(par))-1)^2
}

x0<-c(1,1)
cat("myxp at :")
print(x0)
f0<-myxp(x0, A=aa)
cat("f0 =",f0,"\n")

npar<-2
opxfn<-list2env(list(fn=myxp, gr=NULL, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(A=aa)))

f0u<-ufn(x0,fnuser=opxfn)
cat("f0u =",f0u,"\n")

opxfn$MAXIMIZE=TRUE
f0u<-ufn(x0,fnuser=opxfn)
cat("f0u (MAXIMIZE=TRUE) =",f0u,"\n")

opxfn$FNSCALE=-1
f0u<-ufn(x0,fnuser=opxfn)
cat("f0u (MAXIMIZE=TRUE, FNSCALE=-1) =",f0u,"\n")

opxfn$FNSCALE=0.5
opxfn$MAXIMIZE=FALSE
f0u<-ufn(x0,fnuser=opxfn)
cat("f0u (MAXIMIZE=FALSE, FNSCALE=0.5) =",f0u,"\n")

tmp<-readline("above MAXIMIZE and FNSCALE tests")

axp1<-optim(x0, ufn, control=list(trace=2), fnuser=opxfn)
print(axp1)

cat("=====================================\n\n")

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
npar<-4
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))

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
opxfn<-list2env(list(fn=badlogf, gr=badlogg, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))

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

