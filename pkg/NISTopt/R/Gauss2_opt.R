# Optimization test function Gauss2
# Gauss2 from NISTnls
# ??ref...


Gauss2.f <- function(x) {
   res<-Gauss2.res(x)
   f<-sum(res*res)
}

Gauss2.res <- function(b) {
   xx<-Gauss2$x
   yy<-Gauss2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   res<-b1*exp( -b2*xx ) + b3*exp( -(xx-b4)**2 / b5**2 ) + b6*exp( -(xx-b7)**2 / b8**2 ) - yy
   return(res)
}

# Gauss2 - Jacobian
Gauss2.jac <- function(b) {
   xx<-Gauss2$x
   yy<-Gauss2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   J<-matrix(0,m,n) # define the size of the Jacobian
    expr3 <- exp(-b2 * xx)
    expr5 <- xx - b4
    expr6 <- expr5^2
    expr8 <- b5^2
    expr10 <- exp(-expr6/expr8)
    expr13 <- xx - b7
    expr14 <- expr13^2
    expr16 <- b8^2
    expr18 <- exp(-expr14/expr16)
    J[, 1] <- expr3
    J[, 2] <- -(b1 * (expr3 * xx))
    J[, 3] <- expr10
    J[, 4] <- b3 * (expr10 * (2 * expr5/expr8))
    J[, 5] <- b3 * (expr10 * (expr6 * (2 * b5)/expr8^2))
    J[, 6] <- expr18
    J[, 7] <- b6 * (expr18 * (2 * expr13/expr16))
    J[, 8] <- b6 * (expr18 * (expr14 * (2 * b8)/expr16^2))
   return(J)
}

Gauss2.h <- function(x) {
stop("not defined")
   JJ<-Gauss2.jac(x)
   H <- t(JJ) %*% JJ
   res<-Gauss2.res(x)
stop("not defined")

}

Gauss2.g<-function(x) {
#   stop("not defined")
   JJ<-Gauss2.jac(x)
   res<-Gauss2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Gauss2.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Gauss2) # and load up the data into x and y
   start1 = c( 96,  0.009, 103, 106, 18, 72, 151, 18)
   start2 = c( 98, 0.0105, 103, 105, 20, 73, 150, 20)
}

Gauss2.test<-function() {
}   
