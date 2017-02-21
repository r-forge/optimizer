# Optimization test function misra1d
# misra1d from NISTnls
# ??ref...


misra1d.f <- function(x) {
   res<-misra1d.res(x)
   f<-sum(res*res)
}

misra1d.res <- function(b) {
   xx<-Misra1d$x
   yy<-Misra1d$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1*b2*xx*((1+b2*xx)**(-1)) - yy
   return(res)
}

# misra1d - Jacobian
misra1d.jac <- function(b) {
stop("not defined")
   xx<-misra1d$x
   yy<-misra1d$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

misra1d.h <- function(x) {
stop("not defined")
   JJ<-misra1d.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1d.res(x)
stop("not defined")

}

misra1d.g<-function(x) {
#   stop("not defined")
   JJ<-misra1d.jac(x)
   res<-misra1d.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1d.fgh<-function(x) {
   f<-misra1d.f(x)
   g<-misra1d.g(x)
   H<-misra1d.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

misra1d.setup<-function() {
   library(NISTnls) # get parent collection
   data(Misra1d) # and load up the data into x and y
}

misra1d.test<-function() {
  start1<-c(500, 0.0001)
  start2<-c(450, 0.0003)

}   
