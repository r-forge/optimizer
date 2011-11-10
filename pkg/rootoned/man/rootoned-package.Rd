\name{rootoned-package}
\alias{rootoned-package}
\docType{package}
\title{
rootoned
}
\description{
Roots of functions of one variable.
}
\details{
At 2011-8-18, this package is for finding one real root of a supplied function
of a single variable within a specified interval, or to report that the function
has the same sign at each end of the interval. }
\author{
John C. Nash <nashjc@uottawa.ca>
}
\references{
References to Brent and Dekker need to be added.

Nash J C (1979) Compact numerical methods for computers, Adam-Hilger: Bristol.
Second Edition, 1990, Institute of Physics Publications: Bristol.

}
\keyword{ root-finding }
\keyword{ uniroot }
\seealso{
optimize
}
\examples{
# test rootfinders JN 110817

## Extract to dirctory above rootoned package to test
# source("rootoned/R/root1d.R")
# source("rootoned/R/zeroin.R")

cat("simple test\n")

tf1<-function(x) { 10 - 3*x + x*x*x }
tint<-c(-4,4)
# curve(tf1, tint) ## uncomment to draw

res1<-root1d(tf1, tint)
cat("res1\n")
res1

res2<-zeroin(tf1,tint)
cat("res2\n")
res2


cat("exponential case\n")
alpha<-1.0
efn<-function(x) { exp(-alpha*x) - 0.02 }
tint<-c(0,100)

for (al in 1:10){
   alpha<-0.05*al
#   curve(efn, from=0, to=100)
   resr<-root1d(efn,tint,tol=1e-10)
   cat("root1d:")
   print(resr)
   resz<-zeroin(efn,tint,tol=1e-10)
   cat("zeroin:")
   print(resz)

   tmp<-readline("Next?")
}

cat("Tangent function false root\n")
mytan<-function(xdeg){ # tangent in degrees
    xrad<-xdeg*pi/180.0 # conversion to radians
    tt<-tan(xrad)
}

cat("find root between 80 and 100 degrees\n")
tint<-c(80,100)
rtan<-root1d(mytan, tint)
print(rtan)

}
\keyword{optimize}

