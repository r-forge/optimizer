# J C Nash 2010-2-11 for optimx
# test results in these files are indicated for
# not yet included
# 2011-6-23 -- change from Ravi V. definition to one based on 
#  Nash and Walker-Smith BALF.RES

options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test brown-x.f ...\n")

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

npar<-10 # Down from 500
mres<-npar # number of residuals in sum of squares form
lo<-rep(-100,npar)
hi<-rep(100,npar) # if we use bounds

#p0 <- rnorm(npar,sd=2)
p0=rep(0.5,npar)
system.time(aoptxu <- optimx(par=p0, fn=brown.f, control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]
optansout(aoptxu,filename="./ansbrownu.txt")
aoptxu
system.time(aoptxug <- optimx(par=p0, fn=brown.f, gr=brown.g, control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]
optansout(aoptxug,filename="./ansbrownug.txt")
aoptxug
system.time(aoptxb <- optimx(par=p0, fn=brown.f, control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]
optansout(aoptxb,filename="./ansbrownb.txt")
aoptxb
system.time(aoptxbg <- optimx(par=p0, fn=brown.f, gr=brown.g, control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]
optansout(aoptxbg,filename="./ansbrownbg.txt")
aoptxbg
cat("================= end brown_test ==================\n")
