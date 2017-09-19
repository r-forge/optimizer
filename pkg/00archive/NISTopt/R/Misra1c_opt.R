# Optimization test function Misra1c
# Misra1c from NISTnls
# ??ref...


Misra1c.f <- function(x) {
   res<-Misra1c.res(x)
   f<-sum(res*res)
}

Misra1c.res <- function(b) {
   xx<-Misra1c$x
   yy<-Misra1c$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1 * (1-(1+2*b2*xx)^(-.5)) - yy
   return(res)
}

# Misra1c - Jacobian
Misra1c.jac <- function(b) {
   xx<-Misra1c$x
   yy<-Misra1c$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- 1 + 2 * b2 * xx
   expr4 <- -0.5
   expr6 <- 1 - expr3^expr4
   J[, 1] <- expr6
   J[, 2] <- -(b1 * (expr3^(expr4 - 1) * (expr4 * (2 * xx))))
   return(J)
}

Misra1c.h <- function(x) {
stop("not defined")
   JJ<-Misra1c.jac(x)
   H <- t(JJ) %*% JJ
   res<-Misra1c.res(x)
stop("not defined")

}

Misra1c.g<-function(x) {
#   stop("not defined")
   JJ<-Misra1c.jac(x)
   res<-Misra1c.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Misra1c.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Misra1c) # and load up the data into x and y
   start1<-c(500, 0.0001)
   start2<-c(600, 0.0002)
   out<-list(start1=start1, start2=start2)
   return(out)
}

Misra1c.test<-function() {

}   
