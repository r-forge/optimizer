\name{nls14-package}
\alias{nls14-package}
\alias{nls14}
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
Package: \tab nls14\cr
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

The main functions in \code{nls14} are:

\preformatted{
   nls14fb - Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using a residual and optionally Jacobian
	described as \code{R} functions. 
    20120803: Todo: Make masks more consistent between nls14fb and nls14xb.

   nls14xb - Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using an expression to describe the residual via
        an \code{R} modeling expression. The Jacobian is computed via symbolic
	differentiation.
            
   wrap14nls - Uses nls14xb to solve nonlinear least squares then calls nls() to
            create an object of type nls.

   model2grfun.R - Generate a gradient vector function from a nonlinear 
        model expression and a vector of named parameters.

   model2jacfun.R - Generate a Jacobian matrix function from a nonlinear 
        model expression and a vector of named parameters.

   model2resfun.R - Generate a residual vector function from a nonlinear 
        model expression and a vector of named parameters.

   model2ssfun.R - Generate a sum of squares objective function from a 
        nonlinear model expression and a vector of named parameters.

   modgr.R - compute gradient of the sum of squares function using the 
        Jacobian and residuals for a nonlinear least squares problem
      
   modss.R - computer the sum of squares function from the residuals of
        a nonlinear least squares problem

   myfn.R, mygr.R, myjac.R, myres.R, myss.R - dummy functions that seem to
       be needed so there is an available handle for output of functions that
       generate various functions from expressions.

} 

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
# require(nls14)

traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nls14xb, nls, wrap14nls.
start1  <-  c(b1=1, b2=1, b3=1)
startf1  <-  c(b1=1, b2=1, b3=.1)

eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))

cat("LOCAL DATA IN DATA FRAMES\n")
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weeddata2  <-  data.frame(y=1.5*ydat, tt=tdat)

anls14xb1  <-  try(nls14xb(eunsc, start=start1, trace=traceval, data=weeddata1))
print(anls14xb1)

anls14xb2  <-  try(nls14xb(eunsc, start=start1, trace=traceval, data=weeddata2))
print(anls14xb2)


escal  <-   y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
suneasy  <-  c(b1=200, b2=50, b3=0.3)
ssceasy  <-  c(b1=2, b2=5, b3=3)
st1scal  <-  c(b1=100, b2=10, b3=0.1)



shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
             38.558, 50.156, 62.948, 75.995, 91.972)
    tt  <-  1:12
    res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}
 
shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
    jj  <-  matrix(0.0, 12, 3)
    tt  <-  1:12
    yy  <-  exp(-0.1*x[3]*tt)
    zz  <-  100.0/(1+10.*x[2]*yy)
    jj[tt,1]   <-   zz
    jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
    jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
    return(jj)
}

cat("try nls14fb\n")
st  <-  c(b1=1, b2=1, b3=1)
low  <-  -Inf
up <- Inf


\dontrun{

ans1 <- nls14fb(st, shobbs.res, shobbs.jac, trace=traceval)
ans1
cat("No jacobian function -- use internal approximation\n")
ans1n <- nls14fb(st, shobbs.res, trace=TRUE, control=list(watch=TRUE)) # NO jacfn
ans1n

# tmp <- readline("Try with bounds at 2")
time2 <- system.time(ans2 <- nls14fb(st, shobbs.res, shobbs.jac, upper=c(2,2,2), 
                                  trace=traceval))
ans2
time2


} # end dontrun



cat("BOUNDS")
st2s <- c(b1=1, b2=1, b3=1)

\dontrun{

an1qb1 <- try(nls14xb(escal, start=st2s, trace=traceval, data=weeddata1, 
  lower=c(0,0,0), upper=c(2, 6, 3), control=list(watch=FALSE)))
print(an1qb1)
tmp <- readline("next")
ans2 <- nls14fb(st2s,shobbs.res, shobbs.jac, lower=c(0,0,0), upper=c(2, 6, 3), 
   trace=traceval, control=list(watch=FALSE))
print(ans2)

cat("BUT ... nls() seems to do better from the TRACE information\n")
anlsb <- nls(escal, start=st2s, trace=traceval, data=weeddata1, lower=c(0,0,0),
     upper=c(2,6,3), algorithm='port')
cat("However, let us check the answer\n")
print(anlsb)
cat("BUT...crossprod(resid(anlsb))=",crossprod(resid(anlsb)),"\n")

} # end dontrun


tmp <- readline("next")
cat("Try wrap14nls\n")
traceval <- TRUE
# Data for Hobbs problem
ydat <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat <- seq_along(ydat) # for testing
start1 <- c(b1=1, b2=1, b3=1)
escal <-  y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
up1 <- c(2,6,3)
up2 <- c(1, 5, 9)

weeddata1 <- data.frame(y=ydat, tt=tdat)


\dontrun{

an1w <- try(wrap14nls(escal, start=start1, trace=traceval, data=weeddata1))
print(an1w)



cat("BOUNDED wrap14nls\n")

an1wb <- try(wrap14nls(escal, start=start1, trace=traceval, data=weeddata1, upper=up1))
print(an1wb)


cat("BOUNDED wrap14nls\n")

an2wb <- try(wrap14nls(escal, start=start1, trace=traceval, data=weeddata1, upper=up2))
print(an2wb)

cat("TRY MASKS ONLY\n")

an1xm3 <- try(nls14xb(escal, start1, trace=traceval, data=weeddata1, 
                   masked=c("b3")))
printsum(an1xm3)
an1fm3 <- try(nls14fb(start1, shobbs.res, shobbs.jac, trace=traceval, 
                   data=weeddata1, maskidx=c(3)))
printsum(an1fm3)

an1xm1 <- try(nls14xb
(escal, start1, trace=traceval, data=weeddata1, 
                   masked=c("b1")))
printsum(an1xm1)
an1fm1 <- try(nls14fb(start1, shobbs.res, shobbs.jac, trace=traceval, 
                   data=weeddata1, maskidx=c(1)))
printsum(an1fm1)

} # end dontrun


# Need to check when all parameters masked.??

\dontrun{

cat("\n\n Now check conversion of expression to function\n\n")
cat("K Vandepoel function\n")

x <- c(1,3,5,7) # data
y <- c(37.98,11.68,3.65,3.93)
penetrationks28 <- data.frame(x=x,y=y)

cat("Try nls() -- note the try() function!\n")

fit0  <-  try(nls(y ~ (a+b*exp(1)^(-c * x)), data = penetrationks28, 
    start = c(a=0,b = 1,c=1), trace = TRUE))
print(fit0)

cat("\n\n")

fit1  <-  nls14xb(y ~ (a+b*exp(-c*x)), data = penetrationks28, 
   start = c(a=0,b=1,c=1), trace = TRUE) 
printsum(fit1)

mexprn <- "y ~ (a+b*exp(-c*x))"
pvec <- c(a=0,b=1,c=1)
bnew <- c(a=10,b=3,c=4)

k.r <- model2resfun(mexprn , pvec)
k.j <- model2jacfun(mexprn , pvec)
k.f <- model2ssfun(mexprn , pvec)
k.g <- model2grfun(mexprn , pvec)


cat("At pvec:")
print(pvec)
rp <- k.r(pvec, x=x, y=y)
cat(" rp=")
print(rp)
rf <- k.f(pvec, x=x, y=y)
cat(" rf=")
print(rf)
rj <- k.j(pvec, x=x, y=y)
cat(" rj=")
print(rj)
rg <- k.g(pvec, x=x, y=y)
cat(" rg=")
print(rg)
cat("modss at pvec gives ")
print(modss(pvec, k.r, x=x, y=y))
cat("modgr at pvec gives ")
print(modgr(pvec, k.r, k.j, x=x, y=y))
cat("\n\n")

cat("At bnew:")
print(bnew)
rb <- k.r(bnew, x=x, y=y)
cat(" rb=")
print(rb)
rf <- k.f(bnew, x=x, y=y)
cat(" rf=")
print(rf)
rj <- k.j(bnew, x=x, y=y)
cat(" rj=")
print(rj)
rg <- k.g(bnew, x=x, y=y)
cat(" rg=")
print(rg)
cat("modss at bnew gives ")
print(modss(bnew, k.r, x=x, y=y))
cat("modgr at bnew gives ")
print(modgr(bnew, k.r, k.j, x=x, y=y))
cat("\n\n")

}  ## end of dontrun


}
