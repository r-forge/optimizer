# Optimization test function ENSO
# ENSO from NISTnls
# ??ref...


ENSO.f <- function(x) {
   res<-ENSO.res(x)
   f<-sum(res*res)
}

ENSO.res <- function(b) {
# NOTE: could benefit from some sort of constraint to avoid equal parameters in trig args.
   xx<-ENSO$x # note case!
   yy<-ENSO$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   b9<-b[9]
   res<-b1 + b2*cos( 2*pi*xx/12 ) + b3*sin( 2*pi*xx/12 ) + b5*cos( 2*pi*xx/b4 ) + b6*sin( 2*pi*xx/b4 ) + b8*cos( 2*pi*xx/b7 ) + b9*sin( 2*pi*xx/b7 )  - yy
   return(res)
}

# ENSO - Jacobian
ENSO.jac <- function(b) {
   xx<-ENSO$x
   yy<-ENSO$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   b9<-b[9]
   J<-matrix(0,m,n) # define the size of the Jacobian
   expr2 <- 2 * pi * xx
   expr3 <- expr2/12
   expr4 <- cos(expr3)
   expr7 <- sin(expr3)
   expr10 <- expr2/b4
   expr11 <- cos(expr10)
   expr14 <- sin(expr10)
   expr17 <- expr2/b7
   expr18 <- cos(expr17)
   expr21 <- sin(expr17)
   expr25 <- expr2/b4^2
   expr32 <- expr2/b7^2
   J[, 1] <- 1
   J[, 2] <- expr4
   J[, 3] <- expr7
   J[, 4] <- b5 * (expr14 * expr25) - b6 * (expr11 * expr25)
   J[, 5] <- expr11
   J[, 6] <- expr14
   J[, 7] <- b8 * (expr21 * expr32) - b9 * (expr18 * expr32)
   J[, 8] <- expr18
   J[, 9] <- expr21
   return(J)
}

ENSO.h <- function(x) {
stop("not defined")
   JJ<-ENSO.jac(x)
   H <- t(JJ) %*% JJ
   res<-ENSO.res(x)
stop("not defined")

}

ENSO.g<-function(x) {
   JJ<-ENSO.jac(x)
   res<-ENSO.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

ENSO.setup<-function() {
#   library(NISTnls) # get parent collection
   data(ENSO) # and load up the data into x and y
   start1 = c( 11.0, 3.0, 0.5, 40.0, -0.7, -1.3, 25.0, -0.3, 1.4)
   start2 = c( 10.0, 3.0, 0.5, 44.0, -1.5, 0.5,  26.0, -0.1, 1.5)
   out<-list(start1=start1, start2=start2)
   return(out)
}

ENSO.test<-function() {
}   
