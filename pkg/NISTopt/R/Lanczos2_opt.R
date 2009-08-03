# Optimization test function Lanczos2
# Lanczos2 from NISTnls
# ??ref...


Lanczos2.f <- function(x) {
   res<-Lanczos2.res(x)
   f<-sum(res*res)
}

Lanczos2.res <- function(b) {
   xx<-Lanczos2$x # note case!
   yy<-Lanczos2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   res<-b1*exp(-b2*xx) + b3*exp(-b4*xx) + b5*exp(-b6*xx) - yy
   return(res)
}

# Lanczos2 - Jacobian
Lanczos2.jac <- function(b) {
   xx<-Lanczos2$x
   yy<-Lanczos2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   J<-matrix(0,m,n) # define the size of the Jacobian
    expr3 <- exp(-b2 * xx)
    expr7 <- exp(-b4 * xx)
    expr12 <- exp(-b6 * xx)
    J[, 1] <- expr3
    J[, 2] <- -(b1 * (expr3 * xx))
    J[, 3] <- expr7
    J[, 4] <- -(b3 * (expr7 * xx))
    J[, 5] <- expr12
    J[, 6] <- -(b5 * (expr12 * xx))
    return(J)
}

Lanczos2.h <- function(x) {
stop("not defined")
   JJ<-Lanczos2.jac(x)
   H <- t(JJ) %*% JJ
   res<-Lanczos2.res(x)
stop("not defined")

}

Lanczos2.g<-function(x) {
#   stop("not defined")
   JJ<-Lanczos2.jac(x)
   res<-Lanczos2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Lanczos2.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Lanczos2) # and load up the data into x and y
  start1<-c(1.2,0.3,5.6,5.5,6.5,7.6)
  start2<-c(0.5,0.7,3.6,4.2,4,6.3)
}

Lanczos2.test<-function() {
# ?? fixup
#  start1<-c(1.2,0.3,5.6,5.5,6.5,7.6)
#  start2<-c(0.5,0.7,3.6,4.2,4,6.3)

}   
