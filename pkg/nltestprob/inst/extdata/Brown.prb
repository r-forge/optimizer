# bobbrown.R
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("bobyqa test brown-x.f ...\n")

brown.f <- function(x) {
   res<-brown.res(x)
   fval<-as.numeric(crossprod(res))
}

brown.res <- function(x) { # residuals
   ss<-sum(x)
   nn<-length(x)
   res<-x+(ss-nn-1)*rep(1,nn)
   res[nn]<-prod(x)-1 # Note: do not check for failure!
   res
}

brown.jac <- function(x) { # Jacobian
   nn<-length(x)
   JJ<-matrix(1,nn,nn)+diag(rep(1,nn)) # creates all but last row correctly
   # JJ[nn,]<-1/(x/prod(x))
   JJ[nn,]<-prod(x)/x
   JJ
}

brown.g <- function(x) { # gradient
   res<-brown.res(x) # be nice to keep around and not re-evaluate!
   JJ<-brown.jac(x)
   gg<-as.vector(2*t(JJ)%*%res)
}

source("optimx/R/optimx.R")

npar<-10 # Down from 500
mres<-npar # number of residuals in sum of squares form
lo<-rep(-100,npar)
hi<-rep(100,npar) # if we use bounds

myfn<-list(fn=brown.f, gr=NULL, hess=NULL)


#p0 <- rnorm(npar,sd=2)
p0=rep(0.5,npar)

abo<-optimx(p0, brown.f, lower=lo, upper=hi, method="bobyqa", control=list(trace=2))
##abb<-bobyqa(p0, ufn, lower=lo, upper=hi, control=list(iprint=2), fnuser=myfn)
abb
