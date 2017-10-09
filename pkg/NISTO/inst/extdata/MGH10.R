# Optimization test function mgh10
# mgh10 from NISTnls
# ??ref...


mgh10.f <- function(x) {
   res<-mgh10.res(x)
   f<-sum(res*res)
}

mgh10.res <- function(b) {
   xx<-MGH10$x # note case!
   yy<-MGH10$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<-b1 * exp(b2/(xx+b3)) - yy
   return(res)
}

# mgh10 - Jacobian
mgh10.jac <- function(b) {
stop("not defined")
   xx<-mgh10$x
   yy<-mgh10$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

mgh10.h <- function(x) {
stop("not defined")
   JJ<-mgh10.jac(x)
   H <- t(JJ) %*% JJ
   res<-mgh10.res(x)
stop("not defined")

}

mgh10.g<-function(x) {
#   stop("not defined")
   JJ<-mgh10.jac(x)
   res<-mgh10.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

mgh10.fgh<-function(x) {
   f<-mgh10.f(x)
   g<-mgh10.g(x)
   H<-mgh10.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

mgh10.setup<-function() {
   library(NISTnls) # get parent collection
   data(MGH10) # and load up the data into x and y
}

mgh10.test<-function() {
}   
