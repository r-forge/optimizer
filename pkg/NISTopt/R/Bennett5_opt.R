# Optimization test function Bennett5
# Bennett5 from NISTnls
# ??ref...


Bennett5.f <- function(b) {
   res<-Bennett5.res(b)
   f<-sum(res*res)
}

Bennett5.res <- function(b) {
   xx<-Bennett5$x
   yy<-Bennett5$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-b1*(b2+xx)**(-1/b3) - yy
   return(res)
}

# Bennett5 - Jacobian
Bennett5.jac <- function(b) {
## stop("not defined")
   xx<-Bennett5$x
   yy<-Bennett5$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
#   res<-b1*(b2+xx)**(-1/b3) - yy
    expr1 <- b2 + xx
    expr3 <- -1/b3
    expr4 <- expr1^expr3
    value <- b1 * expr4 - yy
    J[ , 1] <- expr4
    J[ , 2] <- b1 * (expr1^(expr3 - 1) * expr3)
    J[ , 3] <- b1 * (expr4 * (log(expr1) * (1/b3^2)))
    return(J)
}

Bennett5.h <- function(b) {
   JJ<-Bennett5.jac(b)
   H <- t(JJ) %*% JJ
   res<-Bennett5.res(b)
stop("not defined")

}

Bennett5.g<-function(x) {
#   stop("not defined")
   JJ<-Bennett5.jac(x)
   res<-Bennett5.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

Bennett5.setup<-function() {
#   library(NISTnls) # get parent collection
   data(Bennett5) # and load up the data into x and y
   start1 = c( -2000, 50, 0.8)
   start2 = c( -1500, 45, 0.85)
   out<-list(start1=start1, start2=start2)
   return(out)
}

Bennett5.test<-function() {
}   
