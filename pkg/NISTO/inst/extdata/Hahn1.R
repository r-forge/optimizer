# Optimization test function hahn1
# hahn1 from NISTnls
# ??ref...


hahn1.f <- function(x) {
   res<-hahn1.res(x)
   f<-sum(res*res)
}

hahn1.res <- function(b) {
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
   res<-(b1+b2*xx+b3*xx**2+b4*xx**3) / (1+b5*xx+b6*xx**2+b7*xx**3) - yy
   return(res)
}

# hahn1 - Jacobian
hahn1.jac <- function(b) {
stop("not defined")
   xx<-hahn1$x
   yy<-hahn1$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

hahn1.h <- function(x) {
stop("not defined")
   JJ<-hahn1.jac(x)
   H <- t(JJ) %*% JJ
   res<-hahn1.res(x)
stop("not defined")

}

hahn1.g<-function(x) {
#   stop("not defined")
   JJ<-hahn1.jac(x)
   res<-hahn1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

hahn1.fgh<-function(x) {
   f<-hahn1.f(x)
   g<-hahn1.g(x)
   H<-hahn1.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

hahn1.setup<-function() {
   library(NISTnls) # get parent collection
   data(Hahn1) # and load up the data into x and y
}

hahn1.test<-function() {
  start1<-c(10, -1, .05, -0.00001, -0.05, 0.001, -0.000001)
  start2<-c(1, -0.1, .005, -0.000001, -0.005, 0.0001, -0.0000001)


}   
