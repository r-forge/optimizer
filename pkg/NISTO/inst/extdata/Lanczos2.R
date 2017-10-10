# Optimization test function lanczos2
# lanczos2 from NISTnls
# ??ref...


lanczos2.f <- function(x) {
   res<-lanczos2.res(x)
   f<-sum(res*res)
}

lanczos2.res <- function(b) {
   xx<-Lanczos2$x # note case!
   yy<-Lanczos2$y
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

# lanczos2 - Jacobian
lanczos2.jac <- function(b) {
stop("not defined")
   xx<-Lanczos2$x
   yy<-Lanczos2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

lanczos2.h <- function(x) {
stop("not defined")
   JJ<-lanczos2.jac(x)
   H <- t(JJ) %*% JJ
   res<-lanczos2.res(x)
stop("not defined")

}

lanczos2.g<-function(x) {
#   stop("not defined")
   JJ<-lanczos2.jac(x)
   res<-lanczos2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

lanczos2.fgh<-function(x) {
   f<-lanczos2.f(x)
   g<-lanczos2.g(x)
   H<-lanczos2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

lanczos2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Lanczos2) # and load up the data into x and y
}

lanczos2.test<-function() {
  start1<-c(1.2,0.3,5.6,5.5,6.5,7.6)
  start2<-c(0.5,0.7,3.6,4.2,4,6.3)

}   
