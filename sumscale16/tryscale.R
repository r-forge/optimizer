# What is set of parameters so product x maximized for sum(x^2)=1

fx <- function(x) {
    ss <- sqrt(sum(x^2))
    fval<-  - sum( log( x/ss ) ) 
}

fx.g <- function(x){
    nn <- length(x)
    ss <- sqrt(sum(x^2))
    gg <- rep(NA, nn)
    for (ii in 1:nn){ # can vectorize?
       gg[ii] <- 1/x[ii] - 2 * nn* x[ii]/ss   
    }
    gg <- -gg
}

cat("test scaled functions")
x <- 1:6
x <- log(x/10)
cat("enll and gradient\n")
library(numDeriv)
xx <- exp(x)
vx <- fx(xx)
vxg <- fx.g(xx)
vxgn <- grad(fx, xx)
cat("fn, gdiff:", vx, max(abs(vxg-vxgn)),"\n")

