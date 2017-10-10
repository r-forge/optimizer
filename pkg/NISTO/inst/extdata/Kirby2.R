# Optimization test function kirby2
# kirby2 from NISTnls
# ??ref...


kirby2.f <- function(x) {
   res<-kirby2.res(x)
   f<-sum(res*res)
}

kirby2.res <- function(b) {
   xx<-Kirby2$x # note case!
   yy<-Kirby2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   res<-(b1 + b2*xx + b3*xx**2) /(1 + b4*xx + b5*xx**2) - yy
   return(res)
}

# kirby2 - Jacobian
kirby2.jac <- function(b) {
stop("not defined")
   xx<-kirby2$x
   yy<-kirby2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

kirby2.h <- function(x) {
stop("not defined")
   JJ<-kirby2.jac(x)
   H <- t(JJ) %*% JJ
   res<-kirby2.res(x)
stop("not defined")

}

kirby2.g<-function(x) {
#   stop("not defined")
   JJ<-kirby2.jac(x)
   res<-kirby2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

kirby2.fgh<-function(x) {
   f<-kirby2.f(x)
   g<-kirby2.g(x)
   H<-kirby2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

kirby2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Kirby2) # and load up the data into x and y
}

kirby2.test<-function() {
  start1<-c(2,-0.1,0.003, -0.001, 0.00001)
  start2<-c(1.5,-0.15,0.0025, -0.0015, 0.00002)

}   
