\name{nlmrt-package}
\alias{nlmrt-package}
\alias{nlmrt}
\docType{package}
\title{
Tools for solving nonlinear least squares problems.
}
\description{
The package provides some tools related to using the Nash variant
of Levenberg-Marquardt algorithm for nonlinear least squares. 
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
paramters. 

Note: You must provide an expression of the form

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
            
   nlsmnqb - Same approach with bounds and mask constraints.
            
   wrapnls - Uses nlsmnq to solve nonlinear least squares then calls nls() to
            create an object of type nls.
} 

For testing purposes, there are also 
\preformatted{
   nlsmn0 - Nash-Marquardt using solution for delta of 
          (J_transpose*J + lambda*D) * delta = J_transpose*residuals
           via the solve() command.

   nlsmn0b - Same approach with bounds and mask constraints.

   nlsmnc - Nash-Marquardt using solution for delta of 
          (J_transpose*J + lambda*D) * delta = J_transpose*residuals
           via the chol(), backsolve() and forwardsolve() commands.
           We have not (yet!) provided a bounds constrained version of this code.
}

The package contains a number of example or test problems in the \code{inst} 
subdirectory.




}
\author{ John C Nash 

Maintainer:  <nashjc@uottawa.ca>
}
\references{
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     Nash, J. C. and Walker-Smith, M.  (1987) _Nonlinear parameter estimation: 
     an integrated system in BASIC_, Marcel Dekker, New York. Available by
     download from http://macnash.telfer.uottawa.ca/nlpe/

     others!!

}
\keyword{ nls }
\keyword{ nonlinear least squares }
\seealso{\code{nls} }
\examples{
# Want to have:                                        Names of results:
#   - Unconstrained algorithms, inc. nls
#          Simple test with data frame                 anlsmnq1df, 
#                                                      anls1df (singular gradient)
#          Simple test with "global" data              anls1g (singular gradient)
#                                                      anlsmnq1g (must have 'data')
#          -- works with nls, should FAIL with nlmrt
#          Simple test with ... -- NOT INCLUDED
# 
#   - Bounds constrained algorithms -- use data frame approach throughout
#          Repeat Unconstrainted tests
#          Repeat unconstrained nls with port algorithm
#          Simple bounds -- feasible 
#              nlsmn0b 
#              nlsmnqb 
#              nls
#          Repeat with infeasible bounds
#          Masks (only JN codes)
#              1, 2, all parameters
#              special case formulas to check results
#              nls special case formulas
#
#     Put more starts and more problems in test file(s).
#
ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat<-1:length(ydat) # for testing
st1<-c(b1=1, b2=1, b3=1)
eunsc<- y ~ b1/(1+b2*exp(-b3*t))
escal<- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*t))
st1scal<-c(b1=100, b2=10, b3=0.1)
st1easy<-c(b1=200, b2=50, b3=0.3)
st1sceasy<-c(b1=2, b2=5, b3=3)

weeddata<-data.frame(y=ydat, t=tdat)

cat("All 1s start -- unscaled")

# Note that we can put in a table for testing. See inst (rather than tests -- avoid 
#          automatic execution)
anls1df<-try(nls(eunsc, start=st1, trace=TRUE, data=weeddata))
if (class(anls1df) == "try-error") {
   cat("FAILED:")
   print(anls1df)
} else {
   print(anls1df)
}


anlsmnq1df<-nlsmnq(eunsc, start=st1, trace=TRUE, data=weeddata)
print(anlsmnq1df)


y<-ydat
t<-tdat
anls1g<-try(nls(eunsc, start=st1, trace=TRUE))
if (class(anls1g) == "try-error") {
   cat("FAILED:")
   print(anls1g)
} else {
   print(anls1g)
}
anlsmnq1g<-try(nlsmnq(eunsc, start=st1, trace=TRUE))
if (class(anlsmnq1g) == "try-error") {
   cat("FAILED:")
   print(anlsmnq1g)
} else {
   print(anlsmnq1g)
}

\dontrun{
cat("Remove data and try again to see if methods fail\n")
cat("They should complain that data is not available\n")
rm(y)
rm(t)
anls1g<-try(nls(eunsc, start=st1, trace=TRUE))
if (class(anls1g) == "try-error") {
   cat("FAILED:")
   print(anls1g)
} else {
   print(anls1g)
}
anlsmnq1g<-try(nlsmnq(eunsc, start=st1, trace=TRUE))
if (class(anlsmnq1g) == "try-error") {
   cat("FAILED:")
   print(anlsmnq1g)
} else {
   print(anlsmnq1g)
}
} # end dontrun



}
