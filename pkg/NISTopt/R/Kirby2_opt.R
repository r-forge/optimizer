# Optimization test function Kirby2
# Kirby2 from NISTnls
# ??ref...


Kirby2.f <- function(x) {
   res<-Kirby2.res(x)
   f<-sum(res*res)
}

Kirby2.res <- function(b) {
   xx<-Kirby2$x # note case!
   yy<-Kirby2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   res<-(b1 + b2*xx + b3*xx**2) /(1 + b4*xx + b5*xx**2) - yy
   return(res)
}

# Kirby2 - Jacobian
Kirby2.jac <- function(b) {
   xx<-Kirby2$x
   yy<-Kirby2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr3 <- xx*xx
   expr5 <- b1 + b2 * xx + b3 * expr3
   expr9 <- 1 + b4 * xx + b5 * expr3
   expr15 <- expr9*expr9
   J[, 1] <- 1/expr9
   J[, 2] <- xx/expr9
   J[, 3] <- expr3/expr9
   J[, 4] <- -(expr5 * xx/expr15)
   J[, 5] <- -(expr5 * expr3/expr15)
   return(J)
}

Kirby2.h <- function(x) {
stop("not defined")
   JJ<-Kirby2.jac(x)
   H <- t(JJ) %*% JJ
   res<-Kirby2.res(x)
stop("not defined")

}

Kirby2.g<-function(x) {
#   stop("not defined")
   JJ<-Kirby2.jac(x)
   res<-Kirby2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}


Kirby2.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Kirby2) # and load up the data into x and y
  start1<-c(2,-0.1,0.003, -0.001, 0.00001)
  start2<-c(1.5,-0.15,0.0025, -0.0015, 0.00002)
  out<-list(start1=start1, start2=start2)
  return(out)
}

Kirby2.test<-function() {
  # ?? fixup
  start1<-c(2,-0.1,0.003, -0.001, 0.00001)
  start2<-c(1.5,-0.15,0.0025, -0.0015, 0.00002)
}   
