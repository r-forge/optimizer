# Optimization test function Eckerle4
# Eckerle4 from NISTnls
# ??ref...


Eckerle4.f <- function(x) {
   res<-Eckerle4.res(x)
   f<-sum(res*res)
}

Eckerle4.res <- function(b) {
   xx<-Eckerle4$x # note case!
   yy<-Eckerle4$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-(b1/b2) * exp(-0.5*((xx-b3)/b2)**2) - yy
   return(res)
}

# Eckerle4 - Jacobian
Eckerle4.jac <- function(b) {
   xx<-Eckerle4$x
   yy<-Eckerle4$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr1 <- b1/b2
   expr3 <- xx - b3
   expr4 <- expr3/b2
   expr7 <- exp(-0.5 * expr4^2)
   expr9 <- 1/b2
   expr11 <- b2^2
   value <- expr1 * expr7
   J[, 1] <- expr9 * expr7
   J[, 2] <- expr1 * (expr7 * (0.5 * (2 * (expr3/expr11 * expr4)))) - b1/expr11 * expr7
   J[, 3] <- expr1 * (expr7 * (0.5 * (2 * (expr9 * expr4))))
   return(J)
}

Eckerle4.h <- function(x) {
stop("not defined")
   JJ<-Eckerle4.jac(x)
   H <- t(JJ) %*% JJ
   res<-Eckerle4.res(x)
stop("not defined")

}

Eckerle4.g<-function(x) {
#   stop("not defined")
   JJ<-Eckerle4.jac(x)
   res<-Eckerle4.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Eckerle4.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Eckerle4) # and load up the data into x and y
   start1 = c( 1, 10, 500)
   start2 = c( 1.5, 5, 450)
}

Eckerle4.test<-function() {
}   
