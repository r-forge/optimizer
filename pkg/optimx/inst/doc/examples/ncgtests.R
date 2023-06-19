# ncgtests.R
##  author: John C. Nash
rm(list=ls()) # comment out this line if you do not want the workspace cleared
require(optimx)
# source("optimx/R/Rcgminb.R")
sessionInfo()

source("optimx/tests/simplefun.R")

n<-4
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
for (i in 1:n) { 
  lower[i]<-1.0*(i-1)*(n-1)/n
  upper[i]<-1.0*i*(n+1)/n
}
x0<-0.5*(lower+upper)


cat("Now force a mask upper=lower for parameter 3 and see what happens\n")
lower[3] <- upper[3]
x0[3] <- lower[3] # MUST reset parameter also

ncgbdm <- optimr(x0, simfun.f, simfun.g, lower=lower, upper=upper, method="ncg", 
                 control=list(trace=4, watch=TRUE))
proptimr(ncgbdm)

# reset
for (i in 1:n) { 
  lower[i]<-1.0*(i-1)*(n-1)/n
  upper[i]<-1.0*i*(n+1)/n
}
x0<-0.5*(lower+upper)


sncgb <- optimr(x0, fn=simfun.f, gr=simfun.g, lower=lower, upper=upper, method="ncg", 
            control=list(trace=4, maxit=600))
proptimr(sncgb)

sncgqsb <- optimr(x0, fn=simfun.f, gr=simfun.g, lower=lower, upper=upper, method="ncgqs", 
            control=list(trace=4, maxit=600))
proptimr(sncgqsb)

sbvm <- optimr(x0, fn=simfun.f, gr=simfun.g,lower=lower,
               upper=upper, method="Rvmmin", control=list(trace=1))
proptimr(sbvm)


source("optimx/tests/exrosen.R")

n<-4
lo<-rep(-2,n)
up<- -lo # to get arrays set
bdmsk<-rep(1,n)
for (i in 1:n) { 
  lo[i]<-1.2*(i-1)*(n-1)/n
  up[i]<-3*i*(n+1)/n
}
xx<-0.5*(lo+up)
cat("lower:"); print(lo)
cat("upper:"); print(up)
cat("start:"); print(xx)

xncg <- optimr(xx, fn=xrosn.f, gr=xrosn.g, lower=lo, upper=up, method="ncg", control=list(trace=1, maxit=1000))
proptimr(xncg)
xncgqs <- optimr(xx, fn=xrosn.f, gr=xrosn.g, lower=lo, upper=up, method="ncgqs", control=list(trace=1, maxit=1000))
proptimr(xncgqs)

xbvm <- optimr(xx, fn=xrosn.f, gr=xrosn.g,lower=lo,
               upper=up, method="Rvmmin", control=list(trace=1))
proptimr(xbvm)

# xncg0 <- ncg(xx, fn=xrosn.f, gr=xrosn.g, bdmsk=NULL, control=list(trace=4, maxit=600))
# proptimr(xncg0)

xbvm0 <- optimr(xx, fn=xrosn.f, gr=xrosn.g, method="Rvmmin", control=list(trace=1))
proptimr(xbvm0)


xall <- opm(xx, fn=xrosn.f, gr=xrosn.g, lower=lo, upper=up, method="ALL")
summary(xall, order=value)


source("optimx/tests/woodfn.R")
x0 <- c(-3,-1,-3,-1) # Wood standard start
lo <- c(-5, -5, -5, -5)
up <- c(0, 10, 10, 10)

wncg <- ncg(x0, fn=wood.f, gr=wood.g, lower=lo, upper=up, bdmsk=NULL, control=list(trace=4, maxit=600))
proptimr(wncg)

wbvm <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, lower=lo,
               upper=up, method="Rvmmin", control=list(trace=1))
proptimr(wbvm)

