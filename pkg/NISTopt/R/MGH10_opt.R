# Optimization test function MGH10
# MGH10 from NISTnls
# ??ref...


MGH10.f <- function(x) {
   res<-MGH10.res(x)
   f<-sum(res*res)
}

MGH10.res <- function(b) {
   xx<-MGH10$x # note case!
   yy<-MGH10$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<-b1 * exp(b2/(xx+b3)) - yy
   return(res)
}

# MGH10 - Jacobian
MGH10.jac <- function(b) {
   xx<-MGH10$x
   yy<-MGH10$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr1 <- xx + b3
   expr3 <- exp(b2/expr1)
   J[, 1] <- expr3
   J[, 2] <- b1 * (expr3 * (1/expr1))
   J[, 3] <- -(b1 * (expr3 * (b2/(expr1*expr1))))
   return(J)
}

MGH10.h <- function(x) {
stop("not defined")
   JJ<-MGH10.jac(x)
   H <- t(JJ) %*% JJ
   res<-MGH10.res(x)
stop("not defined")

}

MGH10.g<-function(x) {
#   stop("not defined")
   JJ<-MGH10.jac(x)
   res<-MGH10.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

MGH10.setup<-function() {
#   library(NISTnls) # get parent collection
   data(MGH10) # and load up the data into x and y
   start1 = c( 2, 400000, 25000)
   start2 = c( 0.02, 4000, 250)
}

MGH10.test<-function() {
}   
