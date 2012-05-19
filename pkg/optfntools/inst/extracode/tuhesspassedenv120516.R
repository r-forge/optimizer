rm(list=ls())
require(optfntools)
cat("Show how uhess works\n")
source("/home/work/R-optimtest/xdevel/optfntools/R/ugr.R")
source("/home/work/R-optimtest/xdevel/optfntools/R/ufn.R")
source("/home/work/R-optimtest/xdevel/optfntools/R/uhess.R")
## source("/home/work/R-optimtest/develmake/optfntools/R/optstart.R")

# genrosa function code -- matches rosenbrock when npar=2 and gs=100
genrosa.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
        # Note do not at 1.0 so min at 0
	fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

genrosa.g <- function(x, gs=NULL){
# vectorized gradient for genrosa.f
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
#		z2<-1.0 - x[i-1]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

npar<-2
opxfn<-list2env(list(fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0))

x0<-rep(2,2)
hess0<-genrosa.h(x0)

cat("genrose test0\n")
cat("x0:")
print(x0)
cat("hessian from genrosa.h\n")
print(hess0)

cat("hessian from uhess\n")
hess0u<-uhess(x0, opxfn)
print(hess0u)
cat("KFN, KGR, KHESS:",opxfn$KFN, opxfn$KGR, opxfn$KHESS,"\n")
tmp<-readline("Now changed to npar=8")

x1<-rep(pi,8)
npar<-8
opxfn$PARSCALE<-rep(1,8)
hess1<-genrosa.h(x1)

cat("genrose test1\n")
cat("x1:")
print(x1)
cat("hessian from genrosa.h\n")
print(hess1)

cat("hessian from uhess\n")
hess1u<-uhess(x1, opxfn)
print(hess1u)
cat("KFN, KGR, KHESS:",opxfn$KFN, opxfn$KGR, opxfn$KHESS,"\n")
tmp<-readline("Now scale parameters")

opxfn$PARSCALE<-1:npar
cat("PARSCALE:")
print(opxfn$PARSCALE)

x1s<-x1/opxfn$PARSCALE

cat("Scaled x1 = x1s =")
print(x1s)
hess1us<-uhess(x1s, opxfn)
print(hess1us)
rsc<-1/opxfn$PARSCALE
rhess1us<-diag(rsc)%*%hess1us%*%diag(rsc)
cat("Rescaled:\n")
print(rhess1us)
cat("diff=",max(abs(rhess1us-hess1)),"\n")
cat("KFN, KGR, KHESS:",opxfn$KFN, opxfn$KGR, opxfn$KHESS,"\n")

