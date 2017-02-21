# Optimization test function misra1b
# misra1b from NISTnls
# ??ref...


misra1b.f <- function(x) {
   res<-misra1b.res(x)
   f<-sum(res*res)
}

misra1b.res <- function(b) {
   xx<-Misra1b$x
   yy<-Misra1b$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1 * (1-(1+b2*xx/2)**(-2)) - yy
   return(res)
}

# misra1b - Jacobian
misra1b.jac <- function(b) {
stop("not defined")
   xx<-misra1b$x
   yy<-misra1b$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

misra1b.h <- function(x) {
stop("not defined")
   JJ<-misra1b.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1b.res(x)
stop("not defined")

}

misra1b.g<-function(x) {
#   stop("not defined")
   JJ<-misra1b.jac(x)
   res<-misra1b.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1b.fgh<-function(x) {
   f<-misra1b.f(x)
   g<-misra1b.g(x)
   H<-misra1b.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

misra1b.setup<-function() {
   library(NISTnls) # get parent collection
   data(Misra1b) # and load up the data into x and y
}

misra1b.test<-function() {
  start1<-c(500, 0.0001)
  start2<-c(300, 0.0002)

}   
