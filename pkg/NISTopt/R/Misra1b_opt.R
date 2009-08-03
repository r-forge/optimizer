# Optimization test function Misra1b
# Misra1b from NISTnls
# ??ref...


Misra1b.f <- function(x) {
   res<-Misra1b.res(x)
   f<-sum(res*res)
}

Misra1b.res <- function(b) {
   xx<-Misra1b$x
   yy<-Misra1b$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   res<-b1 * (1-(1+b2*xx/2)^(-2)) - yy
   return(res)
}

# Misra1b - Jacobian
Misra1b.jac <- function(b) {
   xx<-Misra1b$x
   yy<-Misra1b$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- 1 + b2 * xx/2
   expr4 <- -2
   expr6 <- 1 - expr3^expr4
   J[, 1] <- expr6
   J[, 2] <- -(b1 * (expr3^(expr4 - 1) * (expr4 * (xx/2))))
   return(J)
}

Misra1b.h <- function(x) {
stop("not defined")
   JJ<-Misra1b.jac(x)
   H <- t(JJ) %*% JJ
   res<-Misra1b.res(x)
stop("not defined")

}

Misra1b.g<-function(x) {
#   stop("not defined")
   JJ<-Misra1b.jac(x)
   res<-Misra1b.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Misra1b.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Misra1b) # and load up the data into x and y
   start1<-c(500, 0.0001)
   start2<-c(300, 0.0002)
}

Misra1b.test<-function() {

}   
