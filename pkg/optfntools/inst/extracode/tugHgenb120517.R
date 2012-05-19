rm(list=ls())
cat("tugHgenb 120517\n")
require(optfntools)

# source("/home/work/R-optimtest/xdevel/optfntools/R/ugHgenb.R")

cat("Rosenbrock, unscaled optimx default\n")

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
trad<-c(-1.2,1)
print(trad)
rf<-fr(trad)
rg<-grr(trad)
print(rf)
print(rg)
npar<-2
opxfn<-list2env(list(fn=fr, gr=grr, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))

# for gs=1 equivalence 20120410
fr1<-function(x){ x1<-x[1]; x2<-x[2]; (x2-x1*x1)^2+(1-x1)^2}

cat("Now the ugHgenb values\n")
ans1<-ugHgenb(trad, fnuser=opxfn, control=list(ktrace=2))
print(ans1)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(rg-ans1$gn)),"\n")
rh<-jacobian(grr, trad)
cat("Hessiant max abs difference: ", max(abs(rh-ans1$Hn)),"\n")
cat("\n\n")
rm(opxfn)
tmp<-readline("now try genrose")

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
#        z2<-1.0 - x[i-1]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
    }
        return(hh)
}

# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
        if(is.null(gs)) { gs=100.0 }
    fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
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
    return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
    n <- length(x)
    hh<-matrix(rep(0, n*n),n,n)
    for (i in 2:n) {
        z1<-x[i]-x[i-1]*x[i-1]
#        z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
    }
        return(hh)
}

trad<-c(-1.2,1)
fval<-genrose.f(trad)
gval<-genrose.g(trad)
Ahess<-genrose.h(trad)
cat("Traditional start\n")
print(trad)
cat("f, g, H\n")
print(fval)
print(gval)
print(Ahess)
cat("\n\n By ufn etc.\n")

myfn<-list2env(list(fn=genrose.f, gr=genrose.g, hess=genrose.h, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))


uf<-ufn(trad, fnuser=myfn)
ugH<-ugHgenb(trad, fnuser=myfn, control=list(ktrace=2))
print(uf)
print(ugH)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gval-ugH$gn)),"\n")
rh<-jacobian(grr, trad)
cat("Hessiant max abs difference: ", max(abs(Ahess-ugH$Hn)),"\n")
cat("\n\n")
rm(myfn)

tmp<-readline("Try alternative genrosa for npar=2 Rosenbrock")
fvala<-genrosa.f(trad)
gvala<-genrosa.g(trad)
Ahessa<-genrosa.h(trad)
cat("Traditional start\n")
print(trad)
npar<-length(trad)
cat("Alt f, g, H\n")
print(fvala)
print(gvala)
print(Ahessa)
cat("\n\n By ufn etc.\n")
myfna<-list2env(list(fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))
ufa<-ufn(trad, fnuser=myfna)
ugHa<-ugHgenb(trad, fnuser=myfna)
print(ufa)
print(ugHa)
gna<-grad(genrosa.f, trad)
hna<-hessian(genrose.f, trad)
rh<-jacobian(grr, trad)
cat("rh:")
print(rh)
cat("numeric grad\n")
print(gna)
cat("numeric hessian\n")
print(hna)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gvala-ugHa$gn)),"\n")
cat("Hessiant max abs difference: ", max(abs(Ahessa-ugHa$Hn)),"\n")
cat("\n\n")
rm(myfna)

tmp<-readline("genrose trad start, but gs=1")
trad<-c(-1.2,1)
fval<-genrosa.f(trad, gs=1)
gval<-genrosa.g(trad, gs=1)
Ahess<-genrosa.h(trad, gs=1)

myfna<-list2env(list(fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(gs=1)))
cat("Traditional start\n")
print(trad)
cat("f, g, H\n")
print(fval)
print(gval)
print(Ahess)
gennog<-ugHgenb(trad,fnuser=myfna)
cat("results of ugHgenb for genrosa at \n")
print(trad)
print(gennog)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gval-gennog$gn)),"\n")
rh<-jacobian(grr, trad)
cat("Hessiant max abs difference: ", max(abs(Ahess-gennog$Hn)),"\n")
cat("\n\n")
rm(myfna)

tmp<-readline("now try higher dimension and different start")

parx<-rep(1,4)
npar<-length(parx)
lower<-rep(-10,4)
upper<-rep(10,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)

myfn<-list2env(list(fn=genrose.f, gr=genrose.g, hess=genrose.h, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))
gennog<-ugHgenb(parx,fnuser=myfn, control=list(ktrace=1))
cat("results of ugHgenb for genrose without gradient code at \n")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gval-gennog$gn)),"\n")
rh<-jacobian(grr, trad)
cat("Hessiant max abs difference: ", max(abs(Ahess-gennog$Hn)),"\n")
cat("*****************************************\n")
cat("\n\n")
rm(myfn)

tmp<-readline("try with hessian set to NULL")

myfn2<-list2env(list(fn=genrose.f, gr=genrose.g, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))
geng<-ugHgenb(parx,fnuser=myfn2, control=list(ktrace=1))
cat("results of ugHgenb for genrose at ")
print(parx)
print(geng)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gval-geng$gn)),"\n")
rh<-jacobian(grr, trad)
cat("Hessiant max abs difference: ", max(abs(Ahess-geng$Hn)),"\n")
cat("*****************************************\n")
cat("\n\n")
rm(myfn2)

tmp<-readline("try from all parameters 0.9, gs=9.4")

parx<-rep(0.9,4)
print(parx)
fval<-genrose.f(parx, gs=9.4)
cat("fn = ",fval,"\n")
gval<-genrose.g(parx, gs=9.4)
cat("g =")
print(gval)
Ahess<-genrose.h(parx, gs=9.4)
cat("Hess =\n")
print(Ahess)

myfn3<-list2env(list(fn=genrose.f, gr=genrose.g, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(gs=9.4)))

gennog<-ugHgenb(parx,fnuser=myfn3, control=list(ktrace=1))

cat("results of ugHgenb with gs=",9.4," for genrose without gradient or Hessian code \n")
print(gennog)
cat("Comparisons\n")
cat("Gradient max abs difference: ", max(abs(gval-gennog$gn)),"\n")
cat("Hessiant max abs difference: ", max(abs(Ahess-gennog$Hn)),"\n")
cat("*****************************************\n")
cat("\n\n")
rm(myfn3)

tmp<-readline("Change gs to 5")
myfn4<-list2env(list(fn=genrose.f, gr=genrose.g, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=list(gs=5)))

cat("\n\nTest with masks and gs=",5,"\n")
msk<-c(1,1,0,1) # masked parameter 3

gengb<-ugHgenb(parx,fnuser=myfn4, bdmsk=msk, control=list(ktrace=1))
print(gengb)
cat("*****************************************\n")

rm(myfn4)

