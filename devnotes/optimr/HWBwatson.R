fnWatson <- function(x) {
        n <- length(x)
        t <- (1:29)/29
        f <- numeric(31)
        for (i in 1:29) {
            f[i] <- sum( (1:(n-1)) * x[2:n] * (t[i]^(0:(n-2))) ) -
                    sum( x * (t[i]^(0:(n-1))) )^2 - 1
        }
        f[30] <- x[1]
        f[31] <- x[2] - x[1]^2 - 1
        sum(f^2)
}

watson.f <- function(x) {
        n <- length(x)
        t <- (1:29)/29
        f <- numeric(31)
        for (i in 1:29) {
            f[i] <- sum( (1:(n-1)) * x[2:n] * (t[i]^(0:(n-2))) ) -
                    sum( x * (t[i]^(0:(n-1))) )^2 - 1
        }
        f[30] <- x[1]
        f[31] <- x[2] - x[1]^2 - 1
        sum(f^2)
}

watson.res <- function(x) {
        n <- length(x)
        t <- (1:29)/29
        f <- numeric(31)
        for (i in 1:29) {
            f[i] <- sum( (1:(n-1)) * x[2:n] * (t[i]^(0:(n-2))) ) -
                    sum( x * (t[i]^(0:(n-1))) )^2 - 1
        }
        f[30] <- x[1]
        f[31] <- x[2] - x[1]^2 - 1
        f
}

library(optimrx)
st6 <- rep(0,6)
## sol6 <- opm(rep(0,6), fnWatson, lower=rep(-2,6), upper=rep(2,6),
##                 method="ALL", control=list(trace=1))
## summary(sol6, order=value)

solfwd6 <- opm(st6, fnWatson, gr="grfwd", lower=rep(-2,6), upper=rep(2,6),
                 method="ALL", control=list(trace=0))
summary(solfwd6, order=value)
  
solfwd6u <- opm(st6, fnWatson, gr="grfwd",
                 method="ALL", control=list(trace=0))
summary(solfwd6u, order=value)

## slbfgs <- optimr(rep(0,6), fnWatson, method="lbfgs", control=list(trace=1))
##  summary(slbfgs)

## ssbplx <- optimr(rep(0,6), fnWatson, method="subplex")
## summary(ssbplx)

library(nlsr)
tnlswat6 <- system.time(nlswat6<-nlfb(st6, watson.res, trace=FALSE))
cat("time nlfb =",tnlswat6,"\n")
summary(nlswat6)

st12 <- rep(0,12)
## sol12 <- opm(rep(0,12), fnWatson, lower=rep(-2,12), upper=rep(2,12),
##                 method="ALL", control=list(trace=1))
## summary(sol12, order=value)

solfwd12 <- opm(st12, fnWatson, gr="grfwd", lower=rep(-2,12), upper=rep(2,12),
                 method="ALL", control=list(trace=0))
summary(solfwd12, order=value)

solfwd12u <- opm(st12, fnWatson, gr="grfwd",
                 method="ALL", control=list(trace=0))
summary(solfwd12u, order=value)

## slbfgs <- optimr(rep(0,12), fnWatson, method="lbfgs", control=list(trace=1))
##  summary(slbfgs)

## ssbplx <- optimr(rep(0,12), fnWatson, method="subplex")
## summary(ssbplx)

library(nlsr)
tnlswat12 <- system.time(nlswat12<-nlfb(st12, watson.res, trace=FALSE))
cat("time nlfb =",tnlswat12,"\n")
summary(nlswat12)
