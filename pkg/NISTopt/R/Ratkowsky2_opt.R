# Optimization test function Ratkowsky2
# Ratkowsky2 from NISTnls
# ??ref...


Ratkowsky2.f <- function(x) {
   res<-Ratkowsky2.res(x)
   f<-sum(res*res)
}

Ratkowsky2.res <- function(b) {
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-b1 / (1+exp(b2-b3*xx)) - yy
   return(res)
}

# Ratkowsky2 - Jacobian
Ratkowsky2.jac <- function(b) {
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- exp(b2 - b3 * xx)
   expr4 <- 1 + expr3
   expr8 <- expr4*expr4
   J[, 1] <- 1/expr4
   J[, 2] <- -(b1 * expr3/expr8)
   J[, 3] <- b1 * (expr3 * xx)/expr8
   return(J)
}

Ratkowsky2.h <- function(x) {
stop("not defined")
   JJ<-Ratkowsky2.jac(x)
   H <- t(JJ) %*% JJ
   res<-Ratkowsky2.res(x)
stop("not defined")

}

Ratkowsky2.g<-function(x) {
#   stop("not defined")
   JJ<-Ratkowsky2.jac(x)
   res<-Ratkowsky2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Ratkowsky2.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Ratkowsky2) # and load up the data into x and y
   start1<-c(100,1,.1)
   start2<-c(75,2.5,0.07)
   start3<-rep(1,3)
   out<-list(start1=start1, start2=start2, start3=start3 )
   return(out)
}

Ratkowsky2.test<-function() {
   start1<-c(100,1,.1)
   start2<-c(75,2.5,0.07)
   start0<-rep(1,3)

}
