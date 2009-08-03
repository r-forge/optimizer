# Optimization test function Nelson
# Nelson from NISTnls
# ??ref...


Nelson.f <- function(x) {
   res<-Nelson.res(x)
   f<-sum(res*res)
}

Nelson.res <- function(b) {
   xx1<-Nelson$x1
   xx2<-Nelson$x2
   yy<-Nelson$y
   logyy<-log(yy)
   res <- rep(NA, length(xx1))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<- b1 - b2*xx1 * exp(-b3*xx2) - logyy
   return(res)
}

# Nelson - Jacobian
Nelson.jac <- function(b) {
   xx1<-Nelson$x1
   xx2<-Nelson$x2
   yy<-Nelson$y
   logyy<-log(yy)
   n<-length(b)
   m<-length(xx1)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr1 <- b2 * xx1
   expr4 <- exp(-b3 * xx2)
   J[, 1] <- 1
   J[, 2] <- -(xx1 * expr4)
   J[, 3] <- expr1 * (expr4 * xx2)
   return(J)
}

Nelson.h <- function(x) {
stop("not defined")
   JJ<-Nelson.jac(x)
   H <- t(JJ) %*% JJ
   res<-Nelson.res(x)
stop("not defined")

}

Nelson.g<-function(x) {
#   stop("not defined")
   JJ<-Nelson.jac(x)
   res<-Nelson.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Nelson.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Nelson) # and load up the data into x and y
   start1<-c(2,0.0001,-0.01)
   start2<-c(2.5, 0.000000005, -0.05)
   start0<-rep(0,3)


}

Nelson.test<-function() {


}   
