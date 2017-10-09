# Optimization test function nelson
# nelson from NISTnls
# ??ref...


nelson.f <- function(x) {
   res<-nelson.res(x)
   f<-sum(res*res)
}

nelson.res <- function(b) {
   xx1<-Nelson$x1
   xx2<-Nelson$x2
   yy<-Nelson$y
   logyy<-log(yy)
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<- b1 - b2*xx1 * exp(-b3*xx2) - logyy
   return(res)
}

# nelson - Jacobian
nelson.jac <- function(b) {
stop("not defined")
   xx1<-Nelson$x1
   xx2<-Nelson$x2
   yy<-Nelson$y
   logyy<-log(yy)
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

nelson.h <- function(x) {
stop("not defined")
   JJ<-nelson.jac(x)
   H <- t(JJ) %*% JJ
   res<-nelson.res(x)
stop("not defined")

}

nelson.g<-function(x) {
#   stop("not defined")
   JJ<-nelson.jac(x)
   res<-nelson.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

nelson.fgh<-function(x) {
   f<-nelson.f(x)
   g<-nelson.g(x)
   H<-nelson.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

nelson.setup<-function() {
   library(NISTnls) # get parent collection
   data(Nelson) # and load up the data into x and y
}

nelson.test<-function() {
  start1<-c(2,0.0001,-0.01)
  start2<-c(2.5, 0.000000005, -0.05)
  start0<-rep(0,3)


}   
