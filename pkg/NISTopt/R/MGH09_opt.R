# Optimization test function MGH09
# MGH09 from NISTnls
# ??ref...


MGH09.f <- function(x) {
   res<-MGH09.res(x)
   f<-sum(res*res)
}

MGH09.res <- function(b) {
   xx<-MGH09$x # note case!
   yy<-MGH09$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<-b1*(xx^2+xx*b2) / (xx^2+xx*b3+b4) - yy
   return(res)
}

# MGH09 - Jacobian
MGH09.jac <- function(b) {
   xx<-MGH09$x
   yy<-MGH09$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   J<-matrix(0,m,n) # define the size of the Jacobian
    expr1 <- xx*xx
    expr3 <- expr1 + xx * b2
    expr4 <- b1 * expr3
    expr7 <- expr1 + xx * b3 + b4
    expr13 <- expr7*expr7
    J[, 1] <- expr3/expr7
    J[, 2] <- b1 * xx/expr7
    J[, 3] <- -(expr4 * xx/expr13)
    J[, 4] <- -(expr4/expr13)
   return(J)
}

MGH09.h <- function(x) {
stop("not defined")
   JJ<-MGH09.jac(x)
   H <- t(JJ) %*% JJ
   res<-MGH09.res(x)
stop("not defined")

}

MGH09.g<-function(x) {
#   stop("not defined")
   JJ<-MGH09.jac(x)
   res<-MGH09.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}


MGH09.setup<-function() {
#   library(NISTnls) # get parent collection
   data(MGH09) # and load up the data into x and y
   start1 = c(25, 39, 41.5, 39)
   start2 = c(0.25, 0.39, 0.415, 0.39)
   out<-list(start1=start1, start2=start2)
   return(out)
}

MGH09.test<-function() {
}   
