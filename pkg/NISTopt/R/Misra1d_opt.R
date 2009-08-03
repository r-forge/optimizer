# Optimization test function Misra1d
# Misra1d from NISTnls
# ??ref...


Misra1d.f <- function(x) {
   res<-Misra1d.res(x)
   f<-sum(res*res)
}

Misra1d.res <- function(b) {
   xx<-Misra1d$x
   yy<-Misra1d$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   res<-b1*b2*xx*((1+b2*xx)**(-1)) - yy
   return(res)
}

# Misra1d - Jacobian
Misra1d.jac <- function(b) {
   xx<-Misra1d$x
   yy<-Misra1d$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr2 <- b1 * b2 * xx
   expr3 <- b2 * xx
   expr4 <- 1 + expr3
   expr5 <- -1
   expr6 <- expr4^expr5
   J[, 1] <- expr3 * expr6
   J[, 2] <- b1 * xx * expr6 + expr2 * (expr4^(expr5 - 1) * (expr5 * xx))
   return(J)
}

Misra1d.h <- function(x) {
stop("not defined")
   JJ<-Misra1d.jac(x)
   H <- t(JJ) %*% JJ
   res<-Misra1d.res(x)
stop("not defined")

}

Misra1d.g<-function(x) {
#   stop("not defined")
   JJ<-Misra1d.jac(x)
   res<-Misra1d.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Misra1d.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Misra1d) # and load up the data into x and y
  start1<-c(500, 0.0001)
  start2<-c(450, 0.0003)
}

Misra1d.test<-function() {

}   
