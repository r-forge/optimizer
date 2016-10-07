## test small resid case with roffset
library(nlsr)

tt <- 1:25
ymod <- 10 * exp(-0.01*tt) + 5
n <- length(tt)
evec0 <- rep(0, n)
evec1 <- 1e-8*runif(n, -.5, .5)
evec2 <-  1e-1*runif(n, -.5, .5)
y0 <- ymod + evec0
y1 <- ymod + evec1
y2 <- ymod + evec2

mydata <- data.frame(tt, y0, y1, y2)

st <- c(aa=1, bb=1, cc=1)

trf <- function(par, data) {
    tt <- data[,"tt"]
    res <- par["aa"]*exp(-par["bb"]*tt) + par["cc"]
}
print(trf(st, data=mydata))
tmp <- readline("cont.")

trj <- function(par, data) {
    tt <- data[,"tt"]
    m <- dim(data)[1]
    JJ <- matrix(NA, nrow=m, ncol=3)
    JJ[,1] <- exp(-par["bb"]*tt)
    JJ[,2] <- - tt * par["aa"] * exp(-par["bb"]*tt)
    JJ[,3] <- 1
    JJ
}
print(trj(st, data=mydata))
tmp <- readline("cont.")

ssf <- function(par, data){
   rr <- trf(par, data)
   ss <- crossprod(rr)
}
library(numDeriv)
print(jacobian(trf, st, data=mydata))
tmp <- readline("cont.")

nlsrfit0 <- try(nlsrxb(y0 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsrfit0
tmp <- readline("cont.")

## following fails in print of nlsrfit00 -- runs out of iterations -- should fix
## nlsrfit00 <- try(nlsrxb(y0 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE, 
##     control=list(rofftest=FALSE, smallsstest=FALSE, lambda=0, laminc=0, lamdec=0)))
## nlsrfit00
tmp <- readline("cont.")

nlsfit0 <-  try(nls(y0 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsfit0
tmp <- readline("cont.")


nlsrfit1 <- try(nlsrxb(y1 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsrfit1
tmp <- readline("cont.")


nlsfit1 <-  try(nls(y1 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsfit1
tmp <- readline("cont.")


nlsrfit2 <- try(nlsrxb(y2 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsrfit2
tmp <- readline("cont.")


nlsfit2 <-  try(nls(y2 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE))
nlsfit2
tmp <- readline("cont.")

fit000<- try(nlsrxb(y0 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE, 
    control=list(rofftest=FALSE, smallsstest=TRUE, lambda=0)))
fit000

fit0000<- try(nlsrxb(y0 ~ aa * exp(-bb*tt) + cc, start=st, data=mydata, trace=TRUE, 
    control=list(watch=TRUE, rofftest=FALSE, smallsstest=TRUE, lambda=0, laminc=0)))
fit0000
