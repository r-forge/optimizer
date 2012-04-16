\name{nlmrt-package}
\alias{nlmrt-package}
\alias{nlmrt}
\docType{package}
\title{
Tools for solving nonlinear least squares problems. UNDER DEVELOPMENT.
}
\description{
The package provides some tools related to using the Nash variant
of Marquardt's algorithm for nonlinear least squares. 
}
\details{
\tabular{ll}{
Package: \tab nlmrt\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-05\cr
License: \tab GPL-2\cr
}
This package includes methods for solving nonlinear least squares problems
specified by a modeling expression and given a starting vector of named 
paramters. Note: You must provide an expression of the form
    lhs ~ rhsexpression
so that the residual expression 
    rhsexpression - lhs
can be computed. The expression can be enclosed in quotes, and this seems to give
fewer difficulties with R. Data variables must already be defined, either within 
the parent environment or else in the dot-arguments. Other symbolic elements in
the modeling expression must be standard functions or else parameters that are 
named in the start vector.

The main functions in \code{nlmrt} are:

\preformatted{
   nlfb - Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using a residual and optionally Jacobian
	described as \code{R} functions. 
    ?? Todo: Make masks more consistent between nlfb and nlxb.

   nlxb - Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using an expression to describe the residual via
        an \code{R} modeling expression. The Jacobian is computed via symbolic
	differentiation.
            
   wrapnls - Uses nlxb to solve nonlinear least squares then calls nls() to
            create an object of type nls.
} 

For testing purposes, there are also some experimental codes using different 
internal computations for the linear algebraic sub-problems in the inst/dev-codes/
sub-folder.
}
\author{ John C Nash 

Maintainer:  <nashjc@uottawa.ca>
}
\references{
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!??

}
\keyword{ nls }
\keyword{ nonlinear least squares }
\seealso{\code{nls} }
\examples{
rm(list=ls())
# library(nlmrt)

# traceval set TRUE to debug or give full history
traceval<-FALSE

# Data for Hobbs problem

ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
y<-ydat  # for testing
tdat<-1:length(ydat) # for testing
# WARNING -- using T can get confusion with TRUE
t<-tdat
# A simple starting vector -- must have named parameters for nlxb, nls, wrapnls.
start1<-c(b1=1, b2=1, b3=1)

eunsc<- y ~ b1/(1+b2*exp(-b3*t))

cat("GLOBAL DATA\n")

anls1g<-try(nls(eunsc, start=start1, trace=traceval))
print(anls1g)

cat("GLOBAL DATA AND EXPRESSION -- SHOULD FAIL\n")
anlxb1g<-try(nlxb(eunsc, start=start1, trace=traceval))
print(anlxb1g)

rm(y)
rm(t)

cat("LOCAL DATA IN DATA FRAMES\n")
weeddata1<-data.frame(y=ydat, t=tdat)
weeddata2<-data.frame(y=1.5*ydat, t=tdat)
anlxb1<-try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1))
print(anlxb1)

anlxb2<-try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata2))
print(anlxb2)

startp1<-c(b1=1, b2=1, b3=.1)

cat("WITH BOUNDS\n")

anlxb1<-try(nlxb(eunsc, start=startp1, lower=c(b1=0, b2=0, b3=0), 
      upper=c(b1=500, b2=100, b3=5), trace=traceval, data=weeddata1))
print(anlxb1)

tmp<-readline("next")

anlsb1<-try(nls(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), 
     upper=c(b1=500, b2=100, b3=5), trace=traceval, data=weeddata1, algorithm='port'))
print(anlsb1)

tmp<-readline("next")

anlxb2<-try(nlxb(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25), trace=traceval, data=weeddata1))
print(anlxb2)


anlsb2<-try(nls(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25), trace=traceval, data=weeddata1, algorithm='port'))
print(anlsb2)

tmp<-readline("next")

cat("TEST MASKS\n")

anlsmnqm<-try(nlxb(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), 
   upper=c(b1=500, b2=100, b3=5), masked=c("b2"), trace=traceval, data=weeddata1))
print(anlsmnqm)

cat("UNCONSTRAINED\n")
an1q<-try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1))
print(an1q)
tmp<-readline("next")


cat("MASKED\n")
an1qm3<-try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1, masked=c("b3")))
print(an1qm3)
tmp<-readline("next")

# Note that the parameters are put in out of order to test code.
an1qm123<-try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1, masked=c("b2","b1","b3")))
print(an1qm123)
tmp<-readline("next")


cat("BOUNDS")
start2<-c(b1=100, b2=10, b3=0.1)
an1qb1<-try(nlxb(eunsc, start=start2, trace=traceval, data=weeddata1, lower=c(0,0,0), upper=c(200, 60, .3)))
print(an1qb1)
tmp<-readline("next")


cat("BOUNDS and MASK")
an1qbm2<-try(nlxb(eunsc, start=start2, trace=traceval, data=weeddata1, lower=c(0,0,0), upper=c(200, 60, .3), masked=c("b2")))
print(an1qbm2)
tmp<-readline("next")


escal<- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*t))
sscaleasy<-c(b1=200, b2=50, b3=0.3)
sunsceasy<-c(b1=2, b2=5, b3=3)
st1scal<-c(b1=100, b2=10, b3=0.1)


cat("EASY start -- unscaled")
anls01<-try(nls(eunsc, start=sscaleasy, trace=traceval, data=weeddata1))
print(anls01)
anlmrt01<-try(nlxb(eunsc, start=sscaleasy, trace=traceval, data=weeddata1))
print(anlmrt01)

cat("All 1s start -- unscaled")
anls02<-try(nls(eunsc, start=start1, trace=traceval, data=weeddata1))
if (class(anls02) == "try-error") {
   cat("FAILED:")
   print(anls02)
} else {
   print(anls02)
}
anlmrt02<-nlxb(eunsc, start=start1, trace=traceval, data=weeddata1)
print(anlmrt02)

cat("ones start -- scaled")
anls03<-try(nls(escal, start=start1, trace=traceval, data=weeddata1))
print(anls03)
anlmrt03<-nlxb(escal, start=start1, trace=traceval, data=weeddata1)
print(anlmrt03)

cat("HARD start -- scaled")
anls04<-try(nls(escal, start=st1scal, trace=traceval, data=weeddata1))
print(anls04)
anlmrt04<-nlxb(escal, start=st1scal, trace=traceval, data=weeddata1)
print(anlmrt04)

cat("EASY start -- scaled")
anls05<-try(nls(escal, start=sscaleasy, trace=traceval, data=weeddata1))
print(anls05)
anlmrt05<-nlxb(escal, start=sscaleasy, trace=traceval, data=weeddata1)
print(anlmrt03)


shobbs.res<-function(x){ # scaled Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    t<-1:12
    res<-100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*t)) - y
}
 
shobbs.jac<-function(x) { # scaled Hobbs weeds problem -- Jacobian
    jj<-matrix(0.0, 12, 3)
    t<-1:12
    yy<-exp(-0.1*x[3]*t)
    zz<-100.0/(1+10.*x[2]*yy)
    jj[t,1] <- zz
    jj[t,2] <- -0.1*x[1]*zz*zz*yy
    jj[t,3] <- 0.01*x[1]*zz*zz*yy*x[2]*t
    return(jj)
}

cat("try nlfb\n")
st<-c(b1=1, b2=1, b3=1)
low<--Inf
up<-Inf

ans1<-nlfb(shobbs.res, shobbs.jac, st, trace=traceval)
ans1

# tmp<-readline("Try with bounds at 2")
up<-c(2,2,2)
time2<-system.time(ans2<-nlfb(shobbs.res, shobbs.jac, st, upper=c(2, 6, 3), trace=traceval))
ans2
time2


cat("BOUNDS")
st2s<-c(b1=1, b2=1, b3=1)
an1qb1<-try(nlxb(escal, start=st2s, trace=traceval, data=weeddata1, 
  lower=c(0,0,0), upper=c(2, 6, 3), control=list(watch=FALSE)))
print(an1qb1)
tmp<-readline("next")
ans2<-nlfb(shobbs.res, shobbs.jac, st2s, lower=c(0,0,0), upper=c(2, 6, 3), 
   trace=traceval, control=list(watch=FALSE))
print(ans2)

cat("BUT ... nls() seems to do better from the TRACE information\n")
anlsb<-nls(escal, start=st2s, trace=traceval, data=weeddata1, lower=c(0,0,0),
     upper=c(2,6,3), algorithm='port')
cat("However, let us check the answer\n")
print(anlsb)
cat("And...crossprod(resid(anlsb))=",crossprod(resid(anlsb)),"\n")

tmp<-readline("next")



an1w<-try(wrapnls(eunsc, start=start1, trace=traceval, data=weeddata1))
print(an1w)

cat("BOUNDED wrapnls\n")

an1wb<-try(wrapnls(eunsc, start=start1, trace=traceval, data=weeddata1, upper=up2))
print(an1wb)


}
