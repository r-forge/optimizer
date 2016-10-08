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

gnjn <- function(start, resfn, jacfn = NULL, trace = FALSE, 
		data=NULL, control=list(), ...){
# simplified Gauss Newton
   offset = 1e6 # for no change in parms
   stepred <- 0.5 # start with this as per nls()
   par <- start
   cat("starting parameters:")
   print(par)
   res <- resfn(par, data, ...)
   ssbest <- as.numeric(crossprod(res))
   cat("initial ss=",ssbest,"\n")
   par0 <- par
   keepon <- TRUE
   while (keepon) {
      cat("SSbest now ", ssbest,"\n")
      JJ <- jacfn(par, data, ...)
      QJ <- qr(JJ)
      delta <- qr.coef(QJ, -res)
      ss <- ssbest + offset*offset # force evaluation
      step <- 1.0
      if (as.numeric(max(par0+delta)+offset) != as.numeric(max(par0+offset)) ) {
         while (ss > ssbest) {
           par <- par0+delta * step
           res <- resfn(par, data, ...)
           ss <- as.numeric(crossprod(res))
           cat("step =", step,"  ss=",ss,"\n")
           print(par)
           tmp <- readline("continue")
           if (ss > ssbest) {
              step <- step * stepred
           } else {
              par0 <- par
              ssbest <- ss
           }
         } # end inner loop
      } else { keepon <- FALSE # done }
   } # end main iteration
} # seems to need this

} # end gnjn

gnjn2 <- function(start, resfn, jacfn = NULL, trace = FALSE, 
		data=NULL, control=list(), ...){
# simplified Gauss Newton
   offset = 1e6 # for no change in parms
   stepred <- 0.5 # start with this as per nls()
   par <- start
   cat("starting parameters:")
   print(par)
   res <- resfn(par, data, ...)
   ssbest <- as.numeric(crossprod(res))
   cat("initial ss=",ssbest,"\n")
   par0 <- par
   keepon <- TRUE
   while (keepon) {
      cat("SSbest now ", ssbest,"\n")
      JJ <- jacfn(par, data, ...)
      JTJ <- crossprod (JJ)
      JTr <- crossprod (JJ, res)
      delta <- - as.vector(solve(JTJ, JTr))
      ss <- ssbest + offset*offset # force evaluation
      step <- 1.0
      if (as.numeric(max(par0+delta)+offset) != as.numeric(max(par0+offset)) ) {
         while (ss > ssbest) {
           par <- par0+delta * step
           res <- resfn(par, data, ...)
           ss <- as.numeric(crossprod(res))
           cat("step =", step,"  ss=",ss,"  best is",ssbest,"\n")
           print(par)
           tmp <- readline("continue")
           if (ss > ssbest) {
              step <- step * stepred
           } else {
              par0 <- par
              ssbest <- ss
           }
         } # end inner loop
      } else { keepon <- FALSE # done }
   } # end main iteration
} # seems to need this
} # end gnjn2


trf <- function(par, data) {
    tt <- data[,"tt"]
    res <- par["aa"]*exp(-par["bb"]*tt) + par["cc"] - y0
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
Ja <- trj(st, data=mydata)
print(Ja)
tmp <- readline("cont.")

library(numDeriv)
Jn <- jacobian(trf, st, data=mydata)
print(Jn)
print(max(abs(Jn-Ja)))
tmp <- readline("cont.")


ssf <- function(par, data){
   rr <- trf(par, data)
   ss <- crossprod(rr)
}

print(ssf(st, data=mydata))
tmp <- readline("cont.")


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
