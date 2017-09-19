# Optimization test function lanczos3
# lanczos3 from NISTnls
# ??ref...


lanczos3.f <- function(x) {
   res<-lanczos3.res(x)
   f<-sum(res*res)
}

lanczos3.res <- function(b) {
   xx<-Lanczos3$x
   yy<-Lanczos3$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   res<-b1*exp(-b2*xx) + b3*exp(-b4*xx) + b5*exp(-b6*xx) - yy
   return(res)
}

# lanczos3 - Jacobian
lanczos3.jac <- function(b) {
stop("not defined")
   xx<-Lanczos3$x
   yy<-Lanczos3$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

lanczos3.h <- function(x) {
stop("not defined")
   JJ<-lanczos3.jac(x)
   H <- t(JJ) %*% JJ
   res<-lanczos3.res(x)
stop("not defined")

}

lanczos3.g<-function(x) {
#   stop("not defined")
   JJ<-lanczos3.jac(x)
   res<-lanczos3.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

lanczos3.fgh<-function(x) {
   f<-lanczos3.f(x)
   g<-lanczos3.g(x)
   H<-lanczos3.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

lanczos3.setup<-function() {
   library(NISTnls) # get parent collection
   data(Lanczos3) # and load up the data into x and y
}

lanczos3.test<-function() {
  start1<-c(1.2,0.3,5.6,5.5,6.5,7.6)
  start2<-c(0.5,0.7,3.6,4.2,4,6.3)

}   
