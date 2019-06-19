pkgname <- "optimr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('optimr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("coef.opm")
### * coef.opm

flush(stderr()); flush(stdout())

### Name: coef.opm
### Title: Summarize opm object
### Aliases: coef.opm coef<-.opm coef<-
### Keywords: nonlinear optimize

### ** Examples

ans <- opm(fn = function(x) sum(x*x), par = 1:2, method="ALL", control=list(trace=1))
coef(ans)

## Not run: 
##D proj <- function(x) x/sum(x)
##D f <- function(x) -prod(proj(x))
##D ans <- opm(1:2, f)
##D ans
##D coef(ans) <- apply(coef(ans), 1, proj)
##D ans
## End(Not run)




cleanEx()
nameEx("hjn")
### * hjn

flush(stderr()); flush(stdout())

### Name: hjn
### Title: Compact R Implementation of Hooke and Jeeves Pattern Search
###   Optimization
### Aliases: hjn
### Keywords: nonlinear optimize

### ** Examples

#####################
## Rosenbrock Banana function
fr <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}

ansrosenbrock0 <- hjn(fn=fr, par=c(1,2), control=list(maxfeval=2000, trace=0))
print(ansrosenbrock0) # use print to allow copy to separate file that 
#    can be called using source()
#####################
# Simple bounds and masks test
bt.f<-function(x){
 sum(x*x)
}

n<-10
xx<-rep(0,n)
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
   lower[i]<-1.0*(i-1)*(n-1)/n
   upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)
ansbt<-hjn(xx, bt.f, lower, upper, bdmsk, control=list(trace=1, maxfeval=2000))

print(ansbt)

#####################
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

xx<-rep(pi,10)
lower<-NULL
upper<-NULL
bdmsk<-NULL
genrosea<-hjn(xx,genrose.f, control=list(maxfeval=2000), gs=10)
print(genrosea)

cat("timings B vs U\n")
lo<-rep(-100,10)
up<-rep(100,10)
bdmsk<-rep(1,10)
tb<-system.time(ab<-hjn(xx,genrose.f, lower=lo, upper=up,
          bdmsk=bdmsk, control=list(trace=0, maxfeval=2000)))[1]
tu<-system.time(au<-hjn(xx,genrose.f, control=list(maxfeval=2000, trace=0)))[1]
cat("times U=",tu,"   B=",tb,"\n")
cat("solution hjnu\n")
print(au)
cat("solution hjnb\n")
print(ab)
cat("diff fu-fb=",au$value-ab$value,"\n")
cat("max abs parameter diff = ", max(abs(au$par-ab$par)),"\n")

maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}




#####################  From Rvmmin.Rd
cat("test bounds and masks\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds1<-hjn(startx,genrose.f, lower=lo,upper=up, control=list(maxfeval=2000, trace=0)) 
print(grbds1)

cat("test lower bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
grbds2<-hjn(startx,genrose.f, lower=lo) 
print(grbds2)

cat("test lower bound single value only\n")
nn<-4
startx<-rep(pi,nn)
lo<-2
up<-rep(10,nn)
grbds3<-hjn(startx,genrose.f, lower=lo) 
print(grbds3)

cat("test upper bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds4<-hjn(startx,genrose.f, upper=up, control=list(maxfeval=2000)) 
print(grbds4)

cat("test upper bound single value only\n")
nn<-4
startx<-rep(pi,nn)
grbds5<-hjn(startx,genrose.f, upper=10, control=list(maxfeval=2000)) 
print(grbds5)



cat("test masks only\n")
nn<-6
bd<-c(1,1,0,0,1,1)
startx<-rep(pi,nn)
grbds6<-hjn(startx,genrose.f, bdmsk=bd, control=list(maxfeval=2000)) 
print(grbds6)

cat("test upper bound on first two elements only\n")
nn<-4
startx<-rep(pi,nn)
upper<-c(10,8, Inf, Inf)
grbds7<-hjn(startx,genrose.f, upper=upper, control=list(maxfeval=2000)) 
print(grbds7)


cat("test lower bound on first two elements only\n")
nn<-4
startx<-rep(0,nn)
lower<-c(0, -0.1 , -Inf, -Inf)
grbds8a<-hjn(startx,genrose.f, lower=lower, control=list(maxfeval=2000)) 
print(grbds8a)

cat("test n=1 problem using simple squares of parameter\n")

sqtst<-function(xx) {
   res<-sum((xx-2)*(xx-2))
}

######### One dimension test
nn<-1
startx<-rep(0,nn)
onepar<-hjn(startx,sqtst,control=list(trace=1)) 
print(onepar)





cleanEx()
nameEx("multistart")
### * multistart

flush(stderr()); flush(stdout())

### Name: multistart
### Title: General-purpose optimization - multiple starts
### Aliases: multistart
### Keywords: nonlinear optimize

### ** Examples

fnR <- function (x, gs=100.0) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(gs * (x1 - x2^2)^2 + (1 - x2)^2)
}
grR <- function (x, gs=100.0) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 4*gs * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 4 * gs * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 2 * gs * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 2 * gs * (x[n] - x[n - 1]^2)
    g
}

pm <- rbind(rep(1,4), rep(pi, 4), rep(-2,4), rep(0,4), rep(20,4))
pm <- as.matrix(pm)
cat("multistart matrix:\n")
print(pm)

ans <- multistart(pm, fnR, grR, method="Rvmmin", control=list(trace=0))
ans




cleanEx()
nameEx("opm")
### * opm

flush(stderr()); flush(stdout())

### Name: opm
### Title: General-purpose optimization
### Aliases: opm
### Keywords: nonlinear optimize

### ** Examples

require(graphics)
cat("Note possible demo(ox) for extended examples\n")


## Show multiple outputs of optimx using all.methods
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

startx<-4*seq(1:10)/3.
ans8<-opm(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
   control=list(all.methods=TRUE, save.failures=TRUE, trace=1), gs=10)
ans8
ans8[, "gevals"]
ans8["spg", ]
summary(ans8, par.select = 1:3)
summary(ans8, order = value)[1, ] # show best value
head(summary(ans8, order = value)) # best few
## head(summary(ans8, order = "value")) # best few -- alternative syntax

## order by value.  Within those values the same to 3 decimals order by fevals.
## summary(ans8, order = list(round(value, 3), fevals), par.select = FALSE)
summary(ans8, order = "list(round(value, 3), fevals)", par.select = FALSE)

## summary(ans8, order = rownames, par.select = FALSE) # order by method name
summary(ans8, order = "rownames", par.select = FALSE) # same

summary(ans8, order = NULL, par.select = FALSE) # use input order
## summary(ans8, par.select = FALSE) # same




cleanEx()
nameEx("optchk")
### * optchk

flush(stderr()); flush(stdout())

### Name: optchk
### Title: General-purpose optimization
### Aliases: optchk
### Keywords: nonlinear optimize

### ** Examples

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

myctrl<- ctrldefault(2)
myctrl$trace <- 3
mychk <- optchk(par=c(-1.2,1), fr, grr, lower=rep(-10,2), upper=rep(10,2), control=myctrl)
cat("result of optchk\n")
print(mychk)




cleanEx()
nameEx("polyopt")
### * polyopt

flush(stderr()); flush(stdout())

### Name: polyopt
### Title: General-purpose optimization - sequential application of methods
### Aliases: polyopt
### Keywords: nonlinear optimize

### ** Examples

fnR <- function (x, gs=100.0) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(gs * (x1 - x2^2)^2 + (1 - x2)^2)
}
grR <- function (x, gs=100.0) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 4*gs * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 4 * gs * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 2 * gs * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 2 * gs * (x[n] - x[n - 1]^2)
    g
}

x0 <- rep(pi, 4)
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(1000, 100), maxfeval= c(1000, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(100, 100), maxfeval= c(100, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans

mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(10, 100), maxfeval= c(10, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans






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
