# Optimization test function MGH17
# MGH17 from NISTnls
# ??ref...


MGH17.f <- function(x) {
   res<-MGH17.res(x)
   f<-sum(res*res)
}

MGH17.res <- function(b) {
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

# MGH17 - Jacobian
MGH17.jac <- function(b) {
   xx<-MGH17$x
   yy<-MGH17$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr1 <- -xx
   expr3 <- exp(expr1 * b4)
   expr7 <- exp(expr1 * b5)
   J[, 1] <- 1
   J[, 2] <- expr3
   J[, 3] <- expr7
   J[, 4] <- -(b2 * (expr3 * xx))
   J[, 5] <- -(b3 * (expr7 * xx))
   return(J)
}

MGH17.h <- function(x) {
stop("not defined")
   JJ<-MGH17.jac(x)
   H <- t(JJ) %*% JJ
   res<-MGH17.res(x)
stop("not defined")

}

MGH17.g<-function(x) {
#   stop("not defined")
   JJ<-MGH17.jac(x)
   res<-MGH17.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}


MGH17.setup<-function() {
#   library(NISTnls) # get parent collection
   data(MGH17) # and load up the data into x and y
   start1 = c( 50,  150, -100, 1,  2)
   start2 = c( 0.5, 1.5, -1,  0.01, 0.02)
   out<-list(start1=start1, start2=start2)
   return(out)
}

MGH17.test<-function() {
}   
