# Optimization test function Roszman1
# Roszman1 from NISTnls
# ??ref...


Roszman1.f <- function(x) {
   res<-Roszman1.res(x)
   f<-sum(res*res)
}

Roszman1.res <- function(b) {
   xx<-Roszman1$x
   yy<-Roszman1$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<- b1 - b2*xx - atan(b3/(xx-b4))/pi  - yy
   return(res)
}

# Roszman1 - Jacobian
Roszman1.jac <- function(b) {
   xx<-Roszman1$x
   yy<-Roszman1$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   # res<- b1 - b2*xx - atan(b3/(xx-b4))/pi  - yy
   J<-matrix(0,m,n) # define the size of the Jacobian
   J[,1]<-1
   J[,2]<--xx
   expr1 <- xx - b4
   expr2 <- b3/expr1
   expr7 <- 1 + expr2^2
   J[,3]<- 1/expr1/expr7/pi
   J[,4]<- b3/expr1^2/expr7/pi
   return(J)
}

Roszman1.h <- function(x) {
stop("not defined")
   JJ<-Roszman1.jac(x)
   H <- t(JJ) %*% JJ
   res<-Roszman1.res(x)
stop("not defined")

}

Roszman1.g<-function(x) {
   JJ<-Roszman1.jac(x)
   res<-Roszman1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Roszman1.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Roszman1) # and load up the data into x and y
  start1<-c(0.1, -0.00001, 1000, -100)
  start2<-c(0.2, -0.000005, 1200, -150)
  start0<-rep(1,4)
}

Roszman1.test<-function() {
  start1<-c(0.1, -0.00001, 1000, -100)
  start2<-c(0.2, -0.000005, 1200, -150)
  start0<-rep(1,4)

NIST<-list()
NIST$value<- 4.9484847331E-04
NIST$par<-c( 2.0196866396E-01,-6.1953516256E-06, 1.2044556708E+03, -1.8134269537E+02)
NIST$ses<-c( 1.9172666023E-02, 3.2058931691E-06, 7.4050983057E+01, 4.9573513849E+01)

}   
