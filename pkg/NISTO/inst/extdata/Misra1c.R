# Optimization test function misra1c
# misra1c from NISTnls
# ??ref...


misra1c.f <- function(x) {
   res<-misra1c.res(x)
   f<-sum(res*res)
}

misra1c.res <- function(b) {
   xx<-Misra1c$x
   yy<-Misra1c$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1 * (1-(1+2*b2*xx)**(-.5)) - yy
   return(res)
}

# misra1c - Jacobian
misra1c.jac <- function(b) {
stop("not defined")
   xx<-misra1c$x
   yy<-misra1c$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

misra1c.h <- function(x) {
stop("not defined")
   JJ<-misra1c.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1c.res(x)
stop("not defined")

}

misra1c.g<-function(x) {
#   stop("not defined")
   JJ<-misra1c.jac(x)
   res<-misra1c.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1c.fgh<-function(x) {
   f<-misra1c.f(x)
   g<-misra1c.g(x)
   H<-misra1c.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

misra1c.setup<-function() {
   library(NISTnls) # get parent collection
   data(Misra1c) # and load up the data into x and y
}

misra1c.test<-function() {
  start1<-c(500, 0.0001)
  start2<-c(600, 0.0002)

}   
