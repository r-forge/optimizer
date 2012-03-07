\name{nlmrt-package}
\alias{nlmrt-package}
\alias{nlmrt}
\docType{package}
\title{
Tools for solving nonlinear least squares problems.
}
\description{
The package provides some tools related to using the Nash variant
of Marquardt's algorithm for nonlinear least squares. 
Todos: different solvers (chol, backsolve, forwardsolve as well
    as qr.solve. Also model2function().
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
   nlsmnq - Nash variant of the Marquardt procedure for nonlinear least squares.
            
   wrapnls - Uses nlsmnq to solve nonlinear least squares then calls nls() to
            create an object of type nls.
} 

For testing purposes, there are also 
\preformatted{
   nlsmn0 - Nash-Marquardt using solution for delta of 
          (J_transpose*J + lambda*D) * delta = J_transpose*residuals
           via the solve() command.

   nlsmnc - Nash-Marquardt using solution for delta of 
          (J_transpose*J + lambda*D) * delta = J_transpose*residuals
           via the chol(), backsolve() and forwardsolve() commands.
}
}
\author{ John C Nash 

Maintainer:  <nashjc@uottawa.ca>
}
\references{
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!

}
\keyword{ nls }
\keyword{ nonlinear least squares }
\seealso{\code{nls} }
\examples{
# to be added

ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat<-1:length(ydat) # for testing
start1<-c(b1=1, b2=1, b3=1)
eunsc<- y ~ b1/(1+b2*exp(-b3*t))
escal<- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*t))
sscaleasy<-c(b1=200, b2=50, b3=0.3)
sunsceasy<-c(b1=2, b2=5, b3=3)
st1scal<-c(b1=100, b2=10, b3=0.1)

weeddata<-data.frame(y=ydat, t=tdat)

cat("EASY start -- unscaled")
anls01<-try(nls(eunsc, start=sscaleasy, trace=TRUE, data=weeddata))
print(anls01)
anlmrt01<-try(nlsmnq(eunsc, start=sscaleasy, trace=TRUE, data=weeddata))
print(anlmrt01)

# ?? can we put in a table of things to run??
cat("All 1s start -- unscaled")
anls02<-try(nls(eunsc, start=start1, trace=TRUE, data=weeddata))
if (class(anls02) == "try-error") {
   cat("FAILED:")
   print(anls02)
} else {
   print(anls02)
}
anlmrt02<-nlsmnq(eunsc, start=start1, trace=TRUE, data=weeddata)
print(anlmrt02)

cat("ones start -- scaled")
anls03<-try(nls(escal, start=start1, trace=TRUE, data=weeddata))
print(anls03)
anlmrt03<-nlsmnq(escal, start=start1, trace=TRUE, data=weeddata)
print(anlmrt03)

cat("HARD start -- scaled")
anls04<-try(nls(escal, start=st1scal, trace=TRUE, data=weeddata))
print(anls04)
anlmrt04<-nlsmnq(escal, start=st1scal, trace=TRUE, data=weeddata)
print(anlmrt04)

cat("EASY start -- scaled")
anls05<-try(nls(escal, start=sscaleasy, trace=TRUE, data=weeddata))
print(anls05)
anlmrt05<-nlsmnq(escal, start=sscaleasy, trace=TRUE, data=weeddata)
print(anlmrt03)



}
