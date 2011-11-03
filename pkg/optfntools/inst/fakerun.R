# fakerun.R
rm(list=ls())
source("fakeopt.R")
#??    require("useroptfn")
fr <- function(x) {   ## Rosenbrock Banana function
#    cat("params:")
#    print(x)
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

userf<-list(fn=fr, gr=grr, hess=NULL)

st1 <- c(-1.2, 1)
lo1 <- c(-5, -5)
up1 <- c(10,10)
bd1 <- c(1,1)
ps1 <- c(1,2)
fs1 <- 2
gs1 <- 100 # NOT USED

require("numDeriv")

cat("using no scaling\n")

f0<-fr(st1)
fu0<-ufn(st1, fnuser=userf, ps=1, fs=1)
cat("f0, fu0:", f0, fu0,"\n")

gua<-ugr(st1, fnuser=userf, ps=1, fs=1)
gaorig<-grr(st1)
gun<-grad(ufn,st1, fnuser=userf, ps=1, fs=1)
gnorig<-grad(fr,st1)
cat("gaorig:")
print(gaorig)
cat("gnorig:")
print(gnorig)
cat("gua:")
print(gua)
cat("gun:")
print(gun)

tmp<-readline("Now with scaling")
f0<-fr(st1)
fu0<-ufn(st1/ps1, fnuser=userf, ps=ps1, fs=fs1)
cat("f0, fu0:", f0, fu0,"\n")

gaorig<-grr(st1)
gnorig<-grad(fr,st1)

gua<-ugr(st1/ps1, fnuser=userf, ps=ps1, fs=fs1)
gun<-grad(ufn,st1/ps1, fnuser=userf, ps=ps1, fs=fs1)

cat("gaorig:")
print(gaorig)
cat("gnorig:")
print(gnorig)
cat("gua:")
print(gua)
cat("gun:")
print(gun)

tmp<-readline("More?")

# need to try different possibilities e.g., gr=NULL here but not in userf
# also no lower or upper etc.
temp<-readline("continue")
#cat("optim with parscale\n")
#aoptf<-optim(st1,fr,method="Nelder-Mead", control=list(parscale=ps1, fnscale=fs1, trace=1))
#aoptf
#temp<-readline("continue")

cat("fakeopt call for a1fg\n")
a1fg <- fakeopt(st1, fn=ufn, gr=ugr, control=list(trace=TRUE), fnuser=userf, 
      ps=ps1, fs=fs1)
print(a1fg)
temp<-readline("continue")

cat("call for a1f\n")
a1f <- fakeopt(st1, fn=ufn, gr=NULL, control=list(trace=TRUE), fnuser=userf, 
      ps=ps1, fs=fs1)
print(a1f)

temp<-readline("continue")

st2 <- rep(pi,4)
lo2 <- rep(-2,4)
up2 <- rep(6,4)
bd2 <- rep(1,4)
ps2 <- seq(1,4)/2
fs2 <- 2
gs2 <- 4



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

userf<-list(fn=genrose.f, gr=genrose.g, hess=genrose.h)

cat("no scaling\n")
f0<-genrose.f(st2, gs=10)
fu0<-ufn(st2, fnuser=userf, ps=1, fs=1, gs=10)
cat("f0, fu0:",f0, fu0,"\n")
cat("gradients:\n")
gaorig<-genrose.g(st2, gs=10)
gnorig<-grad(genrose.f, st2, gs=10)
gua<-ugr(st2, fnuser=userf, ps=1, fs=1, gs=10)
gun<-grad(ufn, st2, fnuser=userf, ps=1, fs=1, gs=10)

cat("gaorig:")
print(gaorig)
cat("gnorig:")
print(gnorig)
cat("gua:")
print(gua)
cat("gun:")
print(gun)

cat("Hessian\n")
haorig<-genrose.h(st2,gs=10)
hnorig<-hessian(genrose.f, st2, gs=10)
hua<-uhess(st2, fnuser=userf, ps=1, fs=1, gs=10)
hun<-hessian(ufn, st2, fnuser=userf, ps=1, fs=1, gs=10)
cat("haorig:")
print(haorig)
cat("hnorig:")
print(hnorig)
cat("Diff = ",max(abs(haorig-hnorig)),"\n")

cat("hua:")
print(hua)
cat("hun:")
print(hun)
cat("Diff = ",max(abs(hua-hun)),"\n")


tmp<-readline("Now with scaling")


cat("With scaling\n")
f0<-genrose.f(st2, gs=10)
fu0<-ufn(st2/ps2, fnuser=userf, ps=ps2, fs=fs2, gs=10)
cat("f0, fu0:",f0, fu0,"\n")
cat("gradients:\n")
gaorig<-genrose.g(st2, gs=10)
gnorig<-grad(genrose.f, st2, gs=10)
gua<-ugr(st2/ps2, fnuser=userf, ps=ps2, fs=fs2, gs=10)
gun<-grad(ufn, st2/ps2, fnuser=userf, ps=ps2, fs=fs2, gs=10)

cat("gaorig:")
print(gaorig)
cat("gnorig:")
print(gnorig)
cat("gua:")
print(gua)
cat("gun:")
print(gun)

cat("Hessian\n")
haorig<-genrose.h(st2,gs=10)
hnorig<-hessian(genrose.f, st2, gs=10)
hua<-uhess(st2/ps2, fnuser=userf, ps=ps2, fs=fs2, gs=10)
hun<-hessian(ufn, st2/ps2, fnuser=userf, ps=ps2, fs=fs2, gs=10)
cat("haorig:")
print(haorig)
cat("hnorig:")
print(hnorig)
cat("Diff = ",max(abs(haorig-hnorig)),"\n")

cat("hua:")
print(hua)
cat("hun:")
print(hun)
cat("Diff = ",max(abs(hua-hun)),"\n")






# cat("call for a2f\n")
# a2fg <- fakeopt(st2, fn=ufn, gr=ugr, control=list(trace=TRUE), fnuser=userf,
#       ps=ps2, fs=fs2, gs=gs2)
# print(a2f)
