# Optimization test function Chwirut2
# Chwirut2 from NISTnls
# ??ref...


Chwirut2.f <- function(x) {
   res<-Chwirut2.res(x)
   f<-sum(res*res)
}

Chwirut2.res <- function(b) {
   xx<-Chwirut2$x # note case !!
   yy<-Chwirut2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-exp(-b1*xx)/(b2+b3*xx) - yy
   return(res)
}

# Chwirut2 - Jacobian
Chwirut2.jac <- function(b) {
## stop("not defined")
   xx<-Chwirut2$x
   yy<-Chwirut2$y
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

Chwirut2.h <- function(x) {
stop("not defined")
   JJ<-Chwirut2.jac(x)
   H <- t(JJ) %*% JJ
   res<-Chwirut2.res(x)
stop("not defined")

}

Chwirut2.g<-function(x) {
#   stop("not defined")
   JJ<-Chwirut2.jac(x)
   res<-Chwirut2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Chwirut2.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Chwirut2) # and load up the data into x and y
# ?? May want to change structure of starts??
   start1 = c(0.1, 0.01, 0.02)
   start2 = c(0.15, 0.008, 0.01)
   out<-list(start1=start1, start2=start2)
   return(out)

}

Chwirut2.test<-function() {
}   
