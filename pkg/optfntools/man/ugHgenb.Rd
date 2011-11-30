% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{ugHgenb}
\alias{ugHgenb}
\encoding{UTF-8}
\title{Generate gradient and Hessian for a function at given parameters.}
\description{
  \code{ugHgenb} is used to generate the gradient and Hessian of an objective
  function used for optimization. If a user-provided gradient function 
  \code{gr} is available it is used to compute the gradient via the wrapper
  \code{ugr}, otherwise 
  package \code{numDeriv} is used. If a user-provided Hessian function
  \code{hess} is available, it is used to compute a Hessian via the wrapper
  \code{uhess}. However, we do not allow the user Hessian function to be
  specified if the user gradient function is NULL.  If the user 
  \code{gr} is available, we use the function \code{jacobian()} from
  package \code{numDeriv} to compute the Hessian. In both these cases we
  check for symmetry of the Hessian. Computational Hessians are commonly
  NOT symmetric. If only the objective function \code{fn} is provided, then
  the Hessian is approximated with the function \code{hessian} from 
  package \code{numDeriv} which guarantees a symmetric matrix. 
}
\usage{
  ugHgenb(par, fnuser, bdmsk=NULL, lower=NULL, upper=NULL,
      control=list(ktrace=0), ps=1, fs=1, maximize=FALSE, numgrad=FALSE, ...) 
}
\arguments{
 \item{par}{Set of parameters, assumed to be at a minimum of the function \code{fn}.}
 \item{fnuser}{Name of the list that has fn=user_objective_function, gr=user_gradient and
        hess=user_hessian. Note that gr or gr and hess may be NULL.}
 \item{bdmsk}{An integer vector of the same length as \code{par}. When an element
     of this vector is 0, the corresponding parameter value is fixed (masked) 
     during an optimization. Non-zero values indicate a parameter is free (1),
     at a lower bound (-3) or at an upper bound (-1), but this routine only
     uses 0 values.}
 \item{lower}{Lower bounds for parameters in \code{par}.}
 \item{upper}{Upper bounds for parameters in \code{par}.}
 \item{control}{A list of controls to the function. Currently 
                  asymptol (default of 1.0e-7 which tests for asymmetry of Hessian approximation
                  (see code for details of the test); 
                  ktrace, a logical flag which, if TRUE, monitors the progress 
                  of ugHgenb (default FALSE), and 
                  stoponerror, defaulting to FALSE to NOT stop when there is
                  an error or asymmetry of Hessian. Set TRUE to stop.}
 \item{ps}{Parameter scaling vector for parameters in \code{par}.}
 \item{fs}{Function scaling value.}
 \item{maximize}{TRUE if we are maximizing a function.}
 \item{numgrad}{TRUE if we are using numerical gradient approximations.}
 \item{\dots}{Extra data needed to compute the function, gradient and Hessian.}
}

\details{
  None
}
\value{
  \code{ansout} a list of four items, 
  \itemize{
      \item{\code{gn}}{  The approximation to the gradient vector.}
      \item{\code{Hn}}{  The approximation to the Hessian matrix.}
      \item{\code{gradOK}}{  TRUE if the gradient has been computed acceptably. FALSE otherwise.}
      \item{\code{hessOK}}{  TRUE if the gradient has been computed acceptably and passes the
             symmetry test. FALSE otherwise.}
      \item{\code{nbm}}{  The number of active bounds and masks.}
  }
}
\examples{

## tugh111029.R -- test ugHgenb.R

require("numDeriv")
#source("/home/john/R-optimtest/work/optfntools/R/ufn.R")
#source("/home/john/R-optimtest/work/optfntools/R/ugr.R")
#source("/home/john/R-optimtest/work/optfntools/R/uhess.R")
#source("/home/john/R-optimtest/work/optfntools/R/ugHgenb.R")

tmp<-readline("Rosenbrock, unscaled optimx default")

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
cat("\n\n")

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
#		z2<-1.0 - x[i-1]
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
#		z2<-1.0-x[i]
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
myfn<-list(fn=genrose.f, gr=genrose.g, hess=genrose.h, callpos=sys.nframe())
uf<-ufn(trad, fnuser=myfn, ps=1, fs=1, maximize=FALSE)
ugH<-ugHgenb(trad, fnuser=myfn, ps=1, fs=1, maximize=FALSE)
print(uf)
print(ugH)


tmp<-readline("Try alternative genrosa that matches Rosenbrock")
fvala<-genrosa.f(trad)
gvala<-genrosa.g(trad)
Ahessa<-genrosa.h(trad)

cat("Traditional start\n")
print(trad)
cat("Alt f, g, H\n")
print(fvala)
print(gvala)
print(Ahessa)
cat("\n\n By ufn etc.\n")
myfn<-list(fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, callpos=sys.nframe())
ufa<-ufn(trad, fnuser=myfn, ps=1, fs=1, maximize=FALSE)
ugHa<-ugHgenb(trad, fnuser=myfn, ps=1, fs=1, maximize=FALSE)
print(ufa)
print(ugHa)
gna<-grad(genrosa.f, trad)
hna<-hessian(genrose.f, trad)
cat("numeric grad\n")
print(gna)
cat("numeric hessian\n")
print(hna)



cat("\n\n")
tmp<-readline("trad start, but gs=1")
trad<-c(-1.2,1)
fval<-genrose.f(trad, gs=1)
gval<-genrose.g(trad, gs=1)
Ahess<-genrose.h(trad, gs=1)
myfn<-list(fn=genrosa.f, gr=genrosa.g, hess=genrosa.h, callpos=sys.nframe())
cat("Traditional start\n")
print(trad)
cat("f, g, H\n")
print(fval)
print(gval)
print(Ahess)
gennog<-ugHgenb(trad,fnuser=myfn, ps=1, fs=1, maximize=FALSE, gs=1)
cat("results of ugHgenb for genrosa at \n")
print(trad)
print(gennog)

cat("\n\n")
tmp<-readline("now try higher dimension and different start")

parx<-rep(1,4)
lower<-rep(-10,4)
upper<-rep(10,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
myfn<-list(fn=genrose.f, gr=genrose.g, hess=genrose.h, callpos=sys.nframe())
gennog<-ugHgenb(parx,fnuser=myfn, control=list(ktrace=1), ps=1, fs=1, maximize=FALSE)
cat("results of ugHgenb for genrose without gradient code at \n")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")

tmp<-readline("try with hessian set to NULL")

myfn2<-list(fn=genrose.f, gr=genrose.g, hess=NULL, callpos=sys.nframe())
geng<-ugHgenb(parx,fnuser=myfn2, control=list(ktrace=1), ps=1, fs=1, maximize=FALSE)
cat("results of ugHgenb for genrose at ")
print(parx)
print(geng)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")

tmp<-readline("try from all parameters 0.9, gs=9.4")

parx<-rep(0.9,4)
fval<-genrose.f(parx, gs=9.4)
gval<-genrose.g(parx, gs=9.4)
Ahess<-genrose.h(parx, gs=9.4)
gennog<-ugHgenb(parx,fnuser=myfn2, control=list(ktrace=1), ps=1, fs=1, maximize=FALSE, gs=9.4)
cat("results of ugHgenb with gs=",9.4," for genrose without gradient code at ")
print(parx)
cat("fn = ",fval,"\n")
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
cat("*****************************************\n")
tmp<-readline("Try call with no ps, fs, maximize (but with fnuser), gs=default")

geng<-ugHgenb(parx, fnuser=myfn2, control=list(ktrace=TRUE))
cat("results of ugHgenb for genrose at ")
print(parx)
print(geng)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)

cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")
tmp<-readline("Change gs to 5")

gst<-5
cat("\n\nTest with masks and gs=",gst,"\n")
msk<-c(1,1,0,1) # masked parameter 3
gengb<-ugHgenb(parx,fnuser=myfn2, bdmsk=msk, control=list(ktrace=TRUE),gs=gst)
print(gengb)

cat("*****************************************\n")


tmp<-readline("load ssfn")

# ===== test scaling of parameters =====

ssfn<-function(par, svec=c(.01, .1, 1, 10, 100), rvec=seq(1,5)){
   sp<-svec*par
   ff<-sum((sp-rvec)^2)
}

##> ff<-ssfn(c(100,20,3,.4,.05))
##> ff #  answer is 0

ssfn.g<-function(par, svec=c(.01, .1, 1, 10, 100), rvec=seq(1,5)){
   sp<-svec*par
   gg<- 2*(sp-rvec)*svec
}

ssfn.h<-function(par, svec=c(.01, .1, 1, 10, 100), rvec=seq(1,5)){
   sp<-svec*par
   hh<-diag(2*svec*svec)
}

tmp<-readline("Try generating ssfn gr and Hessian at solution")

good<-c(100, 20, 3, .4, .05)
ff<-ssfn(good)
cat("At ")
print(good)
cat("ssfn function should be 0 but evaluates to ",ff,"\n")
cat("gradient -- should be small \n")
gg<-ssfn.g(good)
print(gg)
cat("Numerical approximation:\n")
gn<-grad(ssfn,good)
print(gn)
cat("Hessian -- should be diagonal -- evaluates to \n")
hh<-ssfn.h(good)
print(hh)
cat("by numDeriv:\n")
hn<-hessian(ssfn, good)
print(hn)
cat("Max abs diff =", max(abs(hh-hn)),"\n")
tmp<-readline("Now try using jn user functions with scaling")
myfn<-list(fn=ssfn, gr=ssfn.g, hess=ssfn.h, callpos=sys.nframe())
goods<-c(1,2,3,4,5)
psc<-c(100,10,1,.1,.01)
uf<-ufn(goods, fnuser=myfn, fs=1, ps=psc, maximize=FALSE)
cat("uf = ", uf,"\n")
ug<-ugr(goods, fnuser=myfn, fs=1, ps=psc)
cat("ug:")
print(ug)
uh<-uhess(goods, fnuser=myfn, fs=1, ps=psc)
cat("uh\n")
print(uh)
hns<-hessian(ufn, goods, fnuser=myfn, ps=psc)
cat("numerical scaled Hessian:\n")
print(hns)
cat("diff = ",max(abs(uh-hns)),"\n")

tmp<-readline("Now use ugHgenb with good solution")
print(good)
myfn<-list(fn=ssfn, gr=ssfn.g, hess=ssfn.h, callpos=sys.nframe())
ugH0<-ugHgenb(good,fnuser=myfn, control=list(ktrace=1), ps=1, fs=1, maximize=FALSE)
print(ugH0)

tmp<-readline("Now use ugHgenb with goods solution and scaling psc")
print(goods)
print(psc)
ugH1<-ugHgenb(goods,fnuser=myfn, control=list(ktrace=1), ps=psc, fs=1, maximize=FALSE)
print(ugH1)











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




}
\keyword{nonlinear}
\keyword{optimize}
