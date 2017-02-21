# Optimization test function DanielWood
# DanielWood from NISTnls
# ??ref...


DanielWood.f <- function(x) {
   res<-DanielWood.res(x)
   f<-sum(res*res)
}

DanielWood.res <- function(b) {
   xx<-DanielWood$x # case !!
   yy<-DanielWood$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   res<-b1*(xx**b2) - yy
   return(res)
}

# DanielWood - Jacobian
DanielWood.jac <- function(b) {
   xx<-DanielWood$x
   yy<-DanielWood$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr1 <- xx^b2
   J[, 1] <- expr1
   J[, 2] <- b1 * (expr1 * log(xx))
   return(J)
}

DanielWood.h <- function(x) {
stop("not defined")
   JJ<-DanielWood.jac(x)
   H <- t(JJ) %*% JJ
   res<-DanielWood.res(x)
stop("not defined")

}

DanielWood.g<-function(x) {
#   stop("not defined")
   JJ<-DanielWood.jac(x)
   res<-DanielWood.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}


DanielWood.setup<-function() {
#   library(NISTnls) # get parent collection
   data(DanielWood) # and load up the data into x and y
   start1 = c( 1, 5)
   start2 = c( 0.7, 4)
   out<-list(start1=start1, start2=start2)
   return(out)
}

DanielWood.test<-function() {
}   
