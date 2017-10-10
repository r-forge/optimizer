# Optimization test function ratkowsky2
# ratkowsky2 from NISTnls
# ??ref...


ratkowsky2.f <- function(x) {
   res<-ratkowsky2.res(x)
   f<-sum(res*res)
}

ratkowsky2.res <- function(b) {
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-b1 / (1+exp(b2-b3*xx)) - yy
   return(res)
}

# ratkowsky2 - Jacobian
ratkowsky2.jac <- function(b) {
stop("not defined")
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

ratkowsky2.h <- function(x) {
stop("not defined")
   JJ<-ratkowsky2.jac(x)
   H <- t(JJ) %*% JJ
   res<-ratkowsky2.res(x)
stop("not defined")

}

ratkowsky2.g<-function(x) {
#   stop("not defined")
   JJ<-ratkowsky2.jac(x)
   res<-ratkowsky2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

ratkowsky2.fgh<-function(x) {
   f<-ratkowsky2.f(x)
   g<-ratkowsky2.g(x)
   H<-ratkowsky2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

ratkowsky2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Ratkowsky2) # and load up the data into x and y
}

ratkowsky2.test<-function() {
   start1<-c(100,1,.1)
   start2<-c(75,2.5,0.07)
   start0<-rep(1,3)

}
