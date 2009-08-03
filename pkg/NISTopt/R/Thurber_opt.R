# Optimization test function Thurber
# Thurber from NISTnls
# ??ref...


Thurber.f <- function(x) {
   res<-Thurber.res(x)
   f<-sum(res*res)
}

Thurber.res <- function(b) {
   xx<-Thurber$x
   yy<-Thurber$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   res<-(b1 + b2*xx + b3*xx**2 + b4*xx**3) / (1 + b5*xx + b6*xx**2 + b7*xx**3) - yy
   return(res)
}

# Thurber - Jacobian
Thurber.jac <- function(b) {
# stop("not defined")
   xx<-Thurber$x
   yy<-Thurber$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]   
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   J<-matrix(0,m,n) # define the size of the Jacobian
   rn<-(b1 + b2*xx + b3*xx**2 + b4*xx**3) # numerator
   rd<-(1 + b5*xx + b6*xx**2 + b7*xx**3)  # denominator
   fctr<--rn/(rd*rd)
   J[, 1]<-1/rd
   J[, 2]<-xx/rd
   J[, 3]<-J[, 2]*xx
   J[, 4]<-J[, 3]*xx
   J[, 5]<-fctr*xx
   J[, 6]<-J[, 5]*xx
   J[, 7]<-J[, 6]*xx
   return(J)
}

Thurber.h <- function(x) {
stop("not defined")
   JJ<-Thurber.jac(x)
   H <- t(JJ) %*% JJ
   res<-Thurber.res(x)
stop("not defined")

}

Thurber.g<-function(x) {
   stop("not defined")
   JJ<-Thurber.jac(x)
   res<-Thurber.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Thurber.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Thurber) # and load up the data into x and y
}

Thurber.test<-function() {

  start1<-c(1000,1000,400,40,0.7,0.3,0.03)
  start2<-c(1300,1500,500,75,1,0.4,0.05)
  start0<-rep(1,7)

  NIST<-list()
  NIST$value<-5.6427082397E+03
  NIST$par<-c(1.2881396800E+03, 1.4910792535E+03, 5.8323836877E+02, 7.5416644291E+01, 9.6629502864E-01, 3.9797285797E-01,  4.9727297349E-02)
  NIST$ses<-c( 4.6647963344E+00, 3.9571156086E+01,  2.8698696102E+01,  5.5675370270E+00,  3.1333340687E-02, 1.4984928198E-02,  6.5842344623E-03)

}   
