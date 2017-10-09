## @knitr ##PowellSingular.prb
# This is file ##PowellSingular.prb
rm(list=ls())
probname <- "##PowellSingular"
probdesc <- "
Powell Singular function. See More' et al problem 13

Note that the traditional start has integer parameters,
and with a default step of 1, HJ methods do absurdly well.
"

#- Note: environment / list "counters" must already exist

if (exists("pe")) { 
  rm("pe")  
}

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

#- nls format expression
##PowellSingular.formula <- ( y ~ b1*x**b2 )
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
attr(J, "gradient") <- J # needed for nlxb
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

ibrary(numDeriv)

start0 <- c(3, -1, 0, 1) # More et al give f(0,0,0,0) = 0 ?? nonlin equations?
names(start0) <- c("b1","b2","b3","b4")
cat("start0=")
print(start0)
cat("f at start0 =", powsing.f(start0),"\n")
cat("gradient:")
print(powsing.g(start0))
cat("max abs diff from numeric approx =", max(abs(powsing.g(start0)-grad(powsing.f,start0))),"\n")
cat("residuals:")
print(powsing.res(start0))
cat("Jacobian:\n")
print(powsing.jac(start0))
cat("max abs diff from numeric approx =", max(abs(jacobian(powsing.res, start0)-powsing.jac(start0))),"\n")

startpi <- rep(pi,4)
names(startpi) <- c("b1","b2","b3","b4")

library(nlsr)
powsingnlfb0 <- nlfb(start0, powsing.res, powsing.jac, trace=TRUE)
print(powsingnlfb0)

powsingnlfbpi <- nlfb(startpi, powsing.res, powsing.jac, trace=TRUE)
print(powsingnlfbpi)


library(optimr)
powsingopm0 <- opm(start0, powsing.f, powsing.g, method="ALL")
summary(powsingopm0, order=value)

powsingopmpi <- opm(startpi, powsing.f, powsing.g, method="ALL")
summary(powsingopmpi, order=value)

