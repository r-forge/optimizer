# Optimization test function Chwirut1
# Chwirut1 from NISTnls
# ??ref...


Chwirut1.f <- function(x) {
   res<-Chwirut1.res(x)
   f<-sum(res*res)
}

Chwirut1.res <- function(b) {
   xx<-Chwirut1$x # note case!!
   yy<-Chwirut1$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-exp(-b1*xx)/(b2+b3*xx)-yy
   return(res)
}

# Chwirut1 - Jacobian
Chwirut1.jac <- function(b) {
   xx<-Chwirut1$x
   yy<-Chwirut1$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- exp(-b1 * xx)
   expr5 <- b2 + b3 * xx
   expr7 <- expr3 * xx
   expr10 <- expr5*expr5
   value <- expr3/expr5
   J[,1] <- -(expr7/expr5)
   J[,2] <- -(expr3/expr10)
   J[,3] <- -(expr7/expr10)
   return(J)
}

Chwirut1.h <- function(x) {
stop("not defined")
   JJ<-Chwirut1.jac(x)
   H <- t(JJ) %*% JJ
   res<-Chwirut1.res(x)
stop("not defined")

}

Chwirut1.g<-function(x) {
#   stop("not defined")
   JJ<-Chwirut1.jac(x)
   res<-Chwirut1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Chwirut1.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Chwirut1) # and load up the data into x and y
   start1 = c( 0.1, 0.01, 0.02)
   start2 = c( 0.15, 0.008, 0.010)
   out<-list(start1=start1, start2=start2)
   return(out)
}

Chwirut1.test<-function() {
}   
