# Optimization test function Gauss3
# Gauss3 from NISTnls
# ??ref...


Gauss3.f <- function(x) {
   res<-Gauss3.res(x)
   f<-sum(res*res)
}

Gauss3.res <- function(b) {
   xx<-Gauss3$x
   yy<-Gauss3$y
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

# Gauss3 - Jacobian
Gauss3.jac <- function(b) {
stop("not defined")
   xx<-Gauss3$x
   yy<-Gauss3$y
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

Gauss3.h <- function(x) {
stop("not defined")
   JJ<-Gauss3.jac(x)
   H <- t(JJ) %*% JJ
   res<-Gauss3.res(x)
stop("not defined")

}

Gauss3.g<-function(x) {
   JJ<-Gauss3.jac(x)
   res<-Gauss3.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Gauss3.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Gauss3) # and load up the data into x and y
   start1 = c(94.9, 0.009,  90.1, 113, 20, 73.8, 140, 20)
   start2 = c( 96,  0.0096, 80,  110,  25,  74,  139, 25)
}

Gauss3.test<-function() {
  #?? fixup!
  start1<-c(97.0,  0.009,  100.0,  65.0, 20.0,  70.0, 178.0, 16.5)
  start2<-c(94.0, 0.0105,   99.0,  63.0, 25.0,  71.0, 180.0, 20.0)

}   
