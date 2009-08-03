# Optimization test function Ratkowsky3
# Ratkowsky3 from NISTnls
# ??ref...


Ratkowsky3.f <- function(x) {
   res<-Ratkowsky3.res(x)
   f<-sum(res*res)
}

Ratkowsky3.res <- function(b) {
   xx<-Ratkowsky3$x
   yy<-Ratkowsky3$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<- b1 / ((1+exp(b2-b3*xx))**(1/b4)) - yy
   return(res)
}

# Ratkowsky3 - Jacobian
Ratkowsky3.jac <- function(b) {
   xx<-Ratkowsky3$x
   yy<-Ratkowsky3$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   J<-matrix(0,m,n) # define the size of the Jacobian
    expr3 <- exp(b2 - b3 * xx)
    expr4 <- 1 + expr3
    expr5 <- 1/b4
    expr6 <- expr4^expr5
    expr10 <- expr4^(expr5 - 1)
    expr14 <- expr6*expr6
    J[, 1] <- 1/expr6
    J[, 2] <- -(b1 * (expr10 * (expr5 * expr3))/expr14)
    J[, 3] <- b1 * (expr10 * (expr5 * (expr3 * xx)))/expr14
    J[, 4] <- b1 * (expr6 * (log(expr4) * (1/(b4*b4))))/expr14
   return(J)
}

Ratkowsky3.h <- function(x) {
stop("not defined")
   JJ<-Ratkowsky3.jac(x)
   H <- t(JJ) %*% JJ
   res<-Ratkowsky3.res(x)
stop("not defined")

}

Ratkowsky3.g<-function(x) {
#   stop("not defined")
   JJ<-Ratkowsky3.jac(x)
   res<-Ratkowsky3.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Ratkowsky3.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Ratkowsky3) # and load up the data into x and y
   start1<-c(100,10,1,1)
   start2<-c(700,5,0.75,1.3)
   start0<-rep(1,4)
          
}

Ratkowsky3.test<-function() {
   start1<-c(100,10,1,1)
   start2<-c(700,5,0.75,1.3)
   start0<-rep(1,4)
          
NIST<-list()
NIST$value<-8.7864049080E+03
NIST$par<-c(6.9964151270E+02,5.2771253025E+00,7.5962938329E-01,1.2792483859E+00)
NIST$ses<-c(1.6302297817E+01, 2.0828735829E+00, 1.9566123451E-01,6.8761936385E-01)

}
