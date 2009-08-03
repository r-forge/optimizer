# Optimization test function Hahn1
# Hahn1 from NISTnls
# ??ref...


Hahn1.f <- function(x) {
   res<-Hahn1.res(x)
   f<-sum(res*res)
}

Hahn1.res <- function(b) {
   xx<-Hahn1$x # note case!
   yy<-Hahn1$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   res<-(b1+b2*xx+b3*xx^2+b4*xx^3) / (1+b5*xx+b6*xx^2+b7*xx^3) - yy
   return(res)
}

# Hahn1 - Jacobian
Hahn1.jac <- function(b) {
   xx<-Hahn1$x
   yy<-Hahn1$y
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

    expr3 <- xx*xx
    expr6 <- xx*expr3
    expr8 <- b1 + b2 * xx + b3 * expr3 + b4 * expr6
    expr14 <- 1 + b5 * xx + b6 * expr3 + b7 * expr6
    expr21 <- expr14*expr14
    J[, 1] <- 1/expr14
    J[, 2] <- xx/expr14
    J[, 3] <- expr3/expr14
    J[, 4] <- expr6/expr14
    J[, 5] <- -(expr8 * xx/expr21)
    J[, 6] <- -(expr8 * expr3/expr21)
    J[, 7] <- -(expr8 * expr6/expr21)
   return(J)
}

Hahn1.h <- function(x) {
stop("not defined")
   JJ<-Hahn1.jac(x)
   H <- t(JJ) %*% JJ
   res<-Hahn1.res(x)
stop("not defined")

}

Hahn1.g<-function(x) {
#   stop("not defined")
   JJ<-Hahn1.jac(x)
   res<-Hahn1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Hahn1.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Hahn1) # and load up the data into x and y
  start1<-c(10, -1, .05, -0.00001, -0.05, 0.001, -0.000001)
  start2<-c(1, -0.1, .005, -0.000001, -0.005, 0.0001, -0.0000001)
}

Hahn1.test<-function() {
# ?? to fixup
  start1<-c(10, -1, .05, -0.00001, -0.05, 0.001, -0.000001)
  start2<-c(1, -0.1, .005, -0.000001, -0.005, 0.0001, -0.0000001)
}   
