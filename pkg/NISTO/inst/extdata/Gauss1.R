# Optimization test function gauss1
# gauss1 from NISTnls
# ??ref...


gauss1.f <- function(x) {
   res<-gauss1.res(x)
   f<-sum(res*res)
}

gauss1.res <- function(b) {
   xx<-Gauss1$x
   yy<-Gauss1$y
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

# gauss1 - Jacobian
gauss1.jac <- function(b) {
stop("not defined")
   xx<-gauss1$x
   yy<-gauss1$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

gauss1.h <- function(x) {
stop("not defined")
   JJ<-gauss1.jac(x)
   H <- t(JJ) %*% JJ
   res<-gauss1.res(x)
stop("not defined")

}

gauss1.g<-function(x) {
#   stop("not defined")
   JJ<-gauss1.jac(x)
   res<-gauss1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

gauss1.fgh<-function(x) {
   f<-gauss1.f(x)
   g<-gauss1.g(x)
   H<-gauss1.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

gauss1.setup<-function() {
   library(NISTnls) # get parent collection
   data(Gauss1) # and load up the data into x and y
}

gauss1.test<-function() {
  start1<-c(97.0,  0.009,  100.0,  65.0, 20.0,  70.0, 178.0, 16.5)
  start2<-c(94.0, 0.0105,   99.0,  63.0, 25.0,  71.0, 180.0, 20.0)

}   
