pkgname <- "optextras"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('optextras')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("axsearch")
### * axsearch

flush(stderr()); flush(stdout())

### Name: axsearch
### Title: Perform axial search around a supposed minimum and provide
###   diagnostics
### Aliases: axsearch
### Keywords: nonlinear optimize axial search

### ** Examples

#####################
require(optimx)
# Simple bounds test for n=4
bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

n<-4
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
# bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
    lower[i]<-1.0*(i-1)*(n-1)/n
    upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)

cat("lower bounds:")
print(lower)
cat("start:       ")
print(xx)
cat("upper bounds:")
print(upper)

cat("Rvmmin \n\n")
# Note: trace set to 0 below. Change as needed to view progress. 

abtrvm <- optimr(xx, bt.f, bt.g, lower=lower, upper=upper, method="Rvmmin", control=list(trace=0))
# Note: use lower=lower etc. because there is a missing hess= argument
print(abtrvm)

cat("Axial search")
axabtrvm <- axsearch(abtrvm$par, fn=bt.f, fmin=abtrvm$value, lower, upper, bdmsk=NULL, 
              trace=0)
print(axabtrvm)

cat("Now force an early stop\n")
abtrvm1 <- optimr(xx, bt.f, bt.g, lower=lower, upper=upper, method="Rvmmin", 
                  control=list(maxit=1, trace=0))
print(abtrvm1)
cat("Axial search")
axabtrvm1 <- axsearch(abtrvm1$par, fn=bt.f, fmin=abtrvm1$value, lower, upper, bdmsk=NULL, 
                     trace=0)
print(axabtrvm1)


cat("Maximization test\n")
mabtrvm <- optimr(xx, bt.f, bt.g, lower=lower, upper=upper, method="Rvmmin", 
                 control=list(trace=1, maximize=TRUE))
# Note: use lower=lower etc. because there is a missing hess= argument
print(mabtrvm)
cat("Do NOT try axsearch() with maximize\n")
cat("KKT condition check\n")
akktm <- kktchk(mabtrvm$par, bt.f, bt.g, hess=NULL, upper=upper, lower=lower,  maximize=TRUE, control=list(trace=0))
print(akktm)






cleanEx()
nameEx("bmchk")
### * bmchk

flush(stderr()); flush(stdout())

### Name: bmchk
### Title: Check bounds and masks for parameter constraints used in
###   nonlinear optimization
### Aliases: bmchk
### Keywords: nonlinear optimize upper lower bound mask

### ** Examples

#####################

cat("25-dimensional box constrained function\n")
flb <- function(x)
    { p <- length(x); sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) }

start<-rep(2, 25)
cat("\n start:")
print(start)
lo<-rep(2,25)
cat("\n lo:")
print(lo)
hi<-rep(4,25)
cat("\n hi:")
print(hi)
bt<-bmchk(start, lower=lo, upper=hi, trace=1)
print(bt)




cleanEx()
nameEx("bmstep")
### * bmstep

flush(stderr()); flush(stdout())

### Name: bmstep
### Title: Compute the maximum step along a search direction.
### Aliases: bmstep
### Keywords: nonlinear optimize upper lower bound mask

### ** Examples

#####################
xx <- c(1, 1)
lo <- c(0, 0)
up <- c(100, 40)
sdir <- c(4,1)
bm <- c(1,1) # both free
ans <- bmstep(xx, sdir, lo, up, bm, trace=1)
# stepsize
print(ans)
# distance
print(ans*sdir)
# New parameters
print(xx+ans*sdir)




cleanEx()
nameEx("fnchk")
### * fnchk

flush(stderr()); flush(stdout())

### Name: fnchk
### Title: Run tests, where possible, on user objective function
### Aliases: fnchk
### Keywords: optimize

### ** Examples

# Want to illustrate each case.
# Ben Bolker idea for a function that is NOT scalar

benbad<-function(x, y){
   # y may be provided with different structures
   f<-(x-y)^2
} # very simple, but ...

y<-1:10
x<-c(1)
cat("test benbad() with y=1:10, x=c(1)\n")
fc01<-fnchk(x, benbad, trace=1, y)
print(fc01)

y<-as.vector(y)
cat("test benbad() with y=as.vector(1:10), x=c(1)\n")
fc02<-fnchk(x, benbad, trace=1, y)
print(fc02)

y<-as.matrix(y)
cat("test benbad() with y=as.matrix(1:10), x=c(1)\n")
fc03<-fnchk(x, benbad, trace=1, y)
print(fc03)
   
y<-as.array(y)
cat("test benbad() with y=as.array(1:10), x=c(1)\n")
fc04<-fnchk(x, benbad, trace=1, y)
print(fc04)
   
y<-"This is a string"
cat("test benbad() with y a string, x=c(1)\n")
fc05<-fnchk(x, benbad, trace=1, y)
print(fc05)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
xtrad<-c(-1.2,1)
ros1<-fnchk(xtrad, fr, trace=1)
print(ros1)
npar<-2
opros<-list2env(list(fn=fr, gr=NULL, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0, dots=NULL))
uros1<-fnchk(xtrad, ufn, trace=1, fnuser=opros)
print(uros1)
require(optextras)
# require(optimx)
cat("Show how fnchk works\n")




cleanEx()
nameEx("gHgen")
### * gHgen

flush(stderr()); flush(stdout())

### Name: gHgen
### Title: Generate gradient and Hessian for a function at given
###   parameters.
### Aliases: gHgen
### Keywords: nonlinear optimize

### ** Examples

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
#		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

trad<-c(-1.2,1)
ans100fgh<-  gHgen(trad, genrose.f, gr=genrose.g, hess=genrose.h,
      control=list(ktrace=1)) 
print(ans100fgh)
ans100fg<-  gHgen(trad, genrose.f, gr=genrose.g, 
      control=list(ktrace=1)) 
print(ans100fg)
ans100f<-  gHgen(trad, genrose.f, control=list(ktrace=1)) 
print(ans100f)
ans10fgh<-   gHgen(trad, genrose.f, gr=genrose.g, hess=genrose.h,
      control=list(ktrace=1), gs=10) 
print(ans10fgh)
ans10fg<-   gHgen(trad, genrose.f, gr=genrose.g, 
      control=list(ktrace=1), gs=10) 
print(ans10fg)
ans10f<-   gHgen(trad, genrose.f, control=list(ktrace=1), gs=10) 
print(ans10f)




cleanEx()
nameEx("gHgenb")
### * gHgenb

flush(stderr()); flush(stdout())

### Name: gHgenb
### Title: Generate gradient and Hessian for a function at given
###   parameters.
### Aliases: gHgenb
### Keywords: nonlinear optimize

### ** Examples

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
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}


maxfn<-function(x, top=10) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-top-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

parx<-rep(1,4)
lower<-rep(-10,4)
upper<-rep(10,4)
bdmsk<-c(1,1,0,1) # masked parameter 3
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgenb(parx,genrose.f)
cat("results of gHgenb for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgenb(parx,genrose.f,genrose.g)
cat("results of gHgenb for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")
parx<-rep(0.9,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgenb(parx,genrose.f,control=list(ktrace=TRUE), gs=9.4)
cat("results of gHgenb with gs=",9.4," for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgenb(parx,genrose.f,genrose.g, control=list(ktrace=TRUE))
cat("results of gHgenb for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
gst<-5
cat("\n\nTest with full calling sequence and gs=",gst,"\n")
gengall<-gHgenb(parx,genrose.f,genrose.g,genrose.h, control=list(ktrace=TRUE),gs=gst)
print(gengall)


top<-25
x0<-rep(2,4)
cat("\n\nTest for maximization and top=",top,"\n")
cat("Gradient and Hessian will have sign inverted")
maxt<-gHgen(x0, maxfn, control=list(ktrace=TRUE), top=top)
print(maxt)

cat("test against negmaxfn\n")
gneg<-grad(negmaxfn, x0)
Hneg<-hessian(negmaxfn, x0)
# gdiff<-max(abs(gneg-maxt$gn))/max(abs(maxt$gn))
# Hdiff<-max(abs(Hneg-maxt$Hn))/max(abs(maxt$Hn))
# explicitly change sign 
gdiff<-max(abs(gneg-(-1)*maxt$gn))/max(abs(maxt$gn))
Hdiff<-max(abs(Hneg-(-1)*maxt$Hn))/max(abs(maxt$Hn))
cat("gdiff = ",gdiff,"  Hdiff=",Hdiff,"\n")






cleanEx()
nameEx("grback")
### * grback

flush(stderr()); flush(stdout())

### Name: grback
### Title: Backward difference numerical gradient approximation.
### Aliases: grback
### Keywords: optimize

### ** Examples

cat("Example of use of grback\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}

xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grback(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to analytic gradient:\n")
print(ga)

cat("change the step parameter to 1e-4\n")
optsp$deps <- 1e-4
gn2<-grback(xx,myfn, shift=0)
print(gn2)




cleanEx()
nameEx("grcentral")
### * grcentral

flush(stderr()); flush(stdout())

### Name: grcentral
### Title: Central difference numerical gradient approximation.
### Aliases: grcentral
### Keywords: optimize

### ** Examples

cat("Example of use of grcentral\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grcentral(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)




cleanEx()
nameEx("grchk")
### * grchk

flush(stderr()); flush(stdout())

### Name: grchk
### Title: Run tests, where possible, on user objective function and
###   (optionally) gradient and hessian
### Aliases: grchk
### Keywords: optimize

### ** Examples

# Would like examples of success and failure. What about "near misses"??
cat("Show how grchk works\n")
require(optextras)
require(numDeriv)
# require(optimx)

jones<-function(xx){
  x<-xx[1]
  y<-xx[2]
  ff<-sin(x*x/2 - y*y/4)*cos(2*x-exp(y))
  ff<- -ff
}

jonesg <- function(xx) {
  x<-xx[1]
  y<-xx[2]
  gx <-  cos(x * x/2 - y * y/4) * ((x + x)/2) * cos(2 * x - exp(y)) - 
    sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * 2)
  gy <- sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * exp(y)) - cos(x * 
                                                                        x/2 - y * y/4) * ((y + y)/4) * cos(2 * x - exp(y))
  gg <- - c(gx, gy)
}

jonesg2 <- function(xx) {
  gx <- 1
  gy <- 2
  gg <- - c(gx, gy)
}


xx <- c(1, 2)

gcans <- grchk(xx, jones, jonesg, trace=1, testtol=(.Machine$double.eps)^(1/3))
gcans

gcans2 <- grchk(xx, jones, jonesg2, trace=1, testtol=(.Machine$double.eps)^(1/3))
gcans2







cleanEx()
nameEx("grfwd")
### * grfwd

flush(stderr()); flush(stdout())

### Name: grfwd
### Title: Forward difference numerical gradient approximation.
### Aliases: grfwd optsp
### Keywords: optimize

### ** Examples

cat("Example of use of grfwd\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grfwd(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)



cleanEx()
nameEx("grnd")
### * grnd

flush(stderr()); flush(stdout())

### Name: grnd
### Title: A reorganization of the call to numDeriv grad() function.
### Aliases: grnd
### Keywords: nonlinear optimize

### ** Examples

cat("Example of use of grnd\n")
require(numDeriv)
myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grnd(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)



cleanEx()
nameEx("hesschk")
### * hesschk

flush(stderr()); flush(stdout())

### Name: hesschk
### Title: Run tests, where possible, on user objective function and
###   (optionally) gradient and hessian
### Aliases: hesschk
### Keywords: optimize

### ** Examples

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
#		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

trad<-c(-1.2,1)
ans100<-hesschk(trad, genrose.f, genrose.g, genrose.h, trace=1)
print(ans100)
ans10<-hesschk(trad, genrose.f, genrose.g, genrose.h, trace=1, gs=10)
print(ans10)





cleanEx()
nameEx("kktchk")
### * kktchk

flush(stderr()); flush(stdout())

### Name: kktchk
### Title: Check Kuhn Karush Tucker conditions for a supposed function
###   minimum
### Aliases: kktchk
### Keywords: nonlinear optimize

### ** Examples

cat("Show how kktc works\n")

require(optimx)
require(optextras)

jones<-function(xx){
   x<-xx[1]
   y<-xx[2]
   ff<-sin(x*x/2 - y*y/4)*cos(2*x-exp(y))
   ff<- -ff
}

jonesg <- function(xx) {
   x<-xx[1]
   y<-xx[2]
   gx <-  cos(x * x/2 - y * y/4) * ((x + x)/2) * cos(2 * x - exp(y)) - 
         sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * 2)
   gy <- sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * exp(y)) - cos(x * 
          x/2 - y * y/4) * ((y + y)/4) * cos(2 * x - exp(y))
   gg <- - c(gx, gy)
}

xx<-0.5*c(pi,pi)

ans <- optimr(xx, jones, jonesg, method="Rvmmin")
ans

kkans <- kktchk(ans$par, jones, jonesg)
kkans






cleanEx()
nameEx("scalechk")
### * scalechk

flush(stderr()); flush(stdout())

### Name: scalechk
### Title: Check the scale of the initial parameters and bounds input to an
###   optimization code used in nonlinear optimization
### Aliases: scalechk
### Keywords: nonlinear optimize upper lower bound mask

### ** Examples

#####################
  par <- c(-1.2, 1)
  lower <- c(-2, 0)
  upper <- c(100000, 10)
  srat<-scalechk(par, lower, upper,dowarn=TRUE)
  print(srat)
  sratv<-c(srat$lpratio, srat$lbratio)
  if (max(sratv,na.rm=TRUE) > 3) { # scaletol from ctrldefault in optimx
     warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
     cat(warnstr,"\n")
  }




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
