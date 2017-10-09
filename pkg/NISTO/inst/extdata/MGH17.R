# Optimization test function mgh17
# mgh17 from NISTnls
# ??ref...


mgh17.f <- function(x) {
   res<-mgh17.res(x)
   f<-sum(res*res)
}

mgh17.res <- function(b) {
   xx<-MGH17$x # note case!
   yy<-MGH17$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   res<-b1 + b2*exp(-xx*b4) + b3*exp(-xx*b5) - yy
   return(res)
}

# mgh17 - Jacobian
mgh17.jac <- function(b) {
stop("not defined")
   xx<-mgh17$x
   yy<-mgh17$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

mgh17.h <- function(x) {
stop("not defined")
   JJ<-mgh17.jac(x)
   H <- t(JJ) %*% JJ
   res<-mgh17.res(x)
stop("not defined")

}

mgh17.g<-function(x) {
#   stop("not defined")
   JJ<-mgh17.jac(x)
   res<-mgh17.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

mgh17.fgh<-function(x) {
   f<-mgh17.f(x)
   g<-mgh17.g(x)
   H<-mgh17.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

mgh17.setup<-function() {
   library(NISTnls) # get parent collection
   data(MGH17) # and load up the data into x and y
}

mgh17.test<-function() {
}   
