# Optimization test function eckerle4
# eckerle4 from NISTnls
# ??ref...


eckerle4.f <- function(x) {
   res<-eckerle4.res(x)
   f<-sum(res*res)
}

eckerle4.res <- function(b) {
   xx<-Eckerle4$x # note case!
   yy<-Eckerle4$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-(b1/b2) * exp(-0.5*((xx-b3)/b2)**2) - yy
   return(res)
}

# eckerle4 - Jacobian
eckerle4.jac <- function(b) {
stop("not defined")
   xx<-eckerle4$x
   yy<-eckerle4$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

eckerle4.h <- function(x) {
stop("not defined")
   JJ<-eckerle4.jac(x)
   H <- t(JJ) %*% JJ
   res<-eckerle4.res(x)
stop("not defined")

}

eckerle4.g<-function(x) {
#   stop("not defined")
   JJ<-eckerle4.jac(x)
   res<-eckerle4.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

eckerle4.fgh<-function(x) {
   f<-eckerle4.f(x)
   g<-eckerle4.g(x)
   H<-eckerle4.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

eckerle4.setup<-function() {
   library(NISTnls) # get parent collection
   data(Eckerle4) # and load up the data into x and y
}

eckerle4.test<-function() {
}   
