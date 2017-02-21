# Optimization test function mgh09
# mgh09 from NISTnls
# ??ref...


mgh09.f <- function(x) {
   res<-mgh09.res(x)
   f<-sum(res*res)
}

mgh09.res <- function(b) {
   xx<-MGH09$x # note case!
   yy<-MGH09$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<-b1*(xx**2+xx*b2) / (xx**2+xx*b3+b4) - yy
   return(res)
}

# mgh09 - Jacobian
mgh09.jac <- function(b) {
stop("not defined")
   xx<-mgh09$x
   yy<-mgh09$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

mgh09.h <- function(x) {
stop("not defined")
   JJ<-mgh09.jac(x)
   H <- t(JJ) %*% JJ
   res<-mgh09.res(x)
stop("not defined")

}

mgh09.g<-function(x) {
#   stop("not defined")
   JJ<-mgh09.jac(x)
   res<-mgh09.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

mgh09.fgh<-function(x) {
   f<-mgh09.f(x)
   g<-mgh09.g(x)
   H<-mgh09.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

mgh09.setup<-function() {
   library(NISTnls) # get parent collection
   data(MGH09) # and load up the data into x and y
}

mgh09.test<-function() {
}   
