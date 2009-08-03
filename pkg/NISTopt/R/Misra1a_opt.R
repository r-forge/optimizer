# Optimization test function Misra1a
# Misra1a from NISTnls
# ??ref...


Misra1a.f <- function(x) {
   res<-Misra1a.res(x)
   f<-sum(res*res)
}

Misra1a.res <- function(b) {
   xx<-Misra1a$x
   yy<-Misra1a$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1*(1-exp(-b2*xx)) - yy
   return(res)
}

# Misra1a - Jacobian
Misra1a.jac <- function(b) {
   xx<-Misra1a$x
   yy<-Misra1a$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- exp(-b2 * xx)
   expr4 <- 1 - expr3
   J[, 1] <- expr4
   J[, 2] <- b1 * (expr3 * xx)
   return(J)
}

Misra1a.h <- function(x) {
stop("not defined")
   JJ<-Misra1a.jac(x)
   H <- t(JJ) %*% JJ
   res<-Misra1a.res(x)
stop("not defined")

}

Misra1a.g<-function(x) {
#   stop("not defined")
   JJ<-Misra1a.jac(x)
   res<-Misra1a.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Misra1a.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Misra1a) # and load up the data into x and y
  start1<-c(500, 0.0001)
  start2<-c(250, 0.0005)
}

Misra1a.test<-function() {
#  start1<-c(500, 0.0001)
#  start2<-c(250, 0.0005)
}   
