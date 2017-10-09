# Example
# Powell's singular function
#
.Machine$dstep<-1e-7
powsing.res <- function(x) {
rs <- rep(NA, 4)
xx <- x[1]
y <- x[2]
z <- x[3]
w <- x[4]
rs[1] <- xx + 10 * y
rs[2] <- sqrt(5) * (z - w)
rs[3] <- (y - 2*z)^2
rs[4] <- sqrt(10) * (xx - w)^2
rs
}
# Powell's singular function - Jacobian
#
powsing.jac <- function(x) {
xx <- x[1]
y <- x[2]
z <- x[3]
w <- x[4]
J <- matrix(0, 4, 4)
J[1, 1] <-  1
J[1, 2] <-  10
J[2, 3] <-  sqrt(5)
J[2, 4] <-  - sqrt(5)
J[3, 2] <-  2 * (y - 2*z)
J[3, 3] <-  -4 * (y - 2*z)
J[4, 1] <-  2*sqrt(10) * (xx - w)
J[4, 4] <-  -2*sqrt(10) * (xx - w)
J
}

powsing.f<-function(x){
   rs<-powsing.res(x)
   f<-sum(rs^2L)
}

powsing.g<-function(x){
   rs<-powsing.res(x)
   JJ<-powsing.jac(x)
   g<-2*as.numeric(crossprod(JJ,rs))
}

powsing.gf<-function(x){ # NOT EFFICIENT
   f0<-powsing.f(x)
   n<-length(x)
   g<-rep(NA,n)
   ss<-.Machine$dstep
#   cat("dstep=",ss,"\n")
   for (i in 1:n){
       xy<-x
       step<-ss*(abs(x[i])+ss)
       xy[i]<-x[i]+step
       f1<-powsing.f(xy)
       g[i]<-(f1-f0)/step
   }
   g
}   


