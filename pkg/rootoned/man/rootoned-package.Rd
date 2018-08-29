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

Brent

Dekker

Nash J C (1979) Compact numerical methods for computers, Adam-Hilger: Bristol.
Second Edition, 1990, Institute of Physics Publications: Bristol.

}
\keyword{ root-finding }
\keyword{ uniroot }
\seealso{
optimize
}
\examples{
cat("simple polynomial test\n")
tf1<-function(x) { 10 - 3*x + x*x*x }
gtf1<-function(x) { -3 + 3*x*x }
tint<-c(-4,4)
# curve(tf1, tint) ## uncomment to draw

res1<-root1d(tf1, tint)
cat("res1\n")
print(unlist(res1))

res2<-zeroin(tf1,tint)
cat("res2\n")
print(unlist(res2))

x0<-tint[1]
res3a<-newt1d(tf1, gtf1, x0)
print(unlist(res3a))
x0<-tint[2]
res3b<-newt1d(tf1, gtf1, x0)
print(unlist(res3b))



## cat("exponential example\n")

## efn<-function(x) { exp(-alpha*x) - 0.02 }
## gefn<-function(x){-alpha*exp(-alpha*x) }
## mymeth <- c("zeroin", "root1d", "newt1d")
## cat("seems to break here\n")
## for (al in 1:10){
##    alpha<- 0.05*al
##    tint <- c(0, 10/alpha)
##    res<-multrfind(fn=efn, gr=gefn, ri=tint, ftrace=TRUE, meths=mymeth, alpha=alpha)
##    cat("alpha=",alpha,":\n")
##    res
##  }

}
\keyword{optimize}

