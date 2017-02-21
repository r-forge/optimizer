# Optimization test function misra1a
# misra1a from NISTnls
# ??ref...


misra1a.f <- function(x) {
   res<-misra1a.res(x)
   f<-sum(res*res)
}

misra1a.res <- function(b) {
   xx<-Misra1a$x
   yy<-Misra1a$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1*(1-exp(-b2*xx)) - yy
   return(res)
}

# misra1a - Jacobian
misra1a.jac <- function(b) {
stop("not defined")
   xx<-misra1a$x
   yy<-misra1a$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

misra1a.h <- function(x) {
stop("not defined")
   JJ<-misra1a.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1a.res(x)
stop("not defined")

}

misra1a.g<-function(x) {
#   stop("not defined")
   JJ<-misra1a.jac(x)
   res<-misra1a.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1a.fgh<-function(x) {
   f<-misra1a.f(x)
   g<-misra1a.g(x)
   H<-misra1a.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

misra1a.setup<-function() {
   library(NISTnls) # get parent collection
   data(Misra1a) # and load up the data into x and y
}

misra1a.test<-function() {
  start1<-c(500, 0.0001)
  start2<-c(250, 0.0005)

}   
