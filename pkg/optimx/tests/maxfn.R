# maxfn.R
##  author: John C. Nash
rm(list=ls())
cat("maxfn.R -- test that maximize=TRUE works correctly\n")
require(optimx)
sessionInfo()

maxfn<-function(x) {# fn to be MAXIMIZED
  # max = 10 at 1:6
  n<-length(x)
  ss<-seq(1,n)
  f<-10-(crossprod(x-ss))^2
  f<-as.numeric(f)
  return(f)
}

negmaxfn<-function(x) {# explicit negative of maxfn
  f<-(-1)*maxfn(x)
  return(f)
}

n<-4
xx<-rep(1,n) # start at all 1s

maxall <- opm(xx, maxfn, gr="grfwd", method="ALL", control=list(maximize=TRUE))
summary(maxall, order=value)
#?? seems hjkb and hjn fail kkt2
#?? lbfgs does NOT report fevals and gevals

ansmax<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE,trace=2))
print(ansmax) # should work OK

# specifying both maximize and fnscale gives maximize precedence with no message
ansmaxboth<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE, fnscale=-1.0, trace=2))
print(ansmaxboth)

cat("using the negmax function should give same parameters\n")
ansnegmax<-optimr(xx,negmaxfn, gr="grfwd",  method="Rvmmin", control=list(trace=0))
print(ansnegmax)# function value should be -10 however

ansmaxnm<-optimr(xx,maxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=TRUE,trace=2, maxit=2000))
print(ansmaxnm)# try with Nelder-Mead

ansmaxnmbad<-optimr(xx,negmaxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=TRUE,trace=2))
# THIS SHOULD FAIL UNBOUNDED
print(ansmaxnmbad)

# control$maximize takes precedence over control$fnscale=-1.
ansmaxnmgood<-optimr(xx,negmaxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=FALSE,
                            fnscale=-1, trace=2))
print(ansmaxnmgood)


# Test Carlo Lapid suggested fix for optimr()  180313

amaxo<-optimr(xx, maxfn, method="L-BFGS-B", control=list(maximize=TRUE, trace=0))
print(amaxo)

cat("using the negmax function should give same parameters\n")
anegmaxo<-optimr(xx,negmaxfn, method="L-BFGS-B", control=list(trace=0))
print(anegmaxo)



cat("WARNING -- this example should FAIL\n")
cat("maximize=TRUE is NOT set up in hjn()\n")
# 160706 -- not set up to maximize, except through optimr perhaps
n<-4
xx<-rep(1,n)
ansmax<-try(hjn(xx,maxfn, control=list(maximize=TRUE,trace=1, maxit=10, maxfeval=2000)))
print(ansmax)

cat("\nTry to maximize through optimr()\n")
ansmaxo<-try(optimr(xx,maxfn, method="hjn", control=list(maximize=TRUE,trace=1, maxit=10, maxfeval=2000)))
print(ansmaxo)

cat("using the negmax function should give same parameters\n")
ansnegmax<-hjn(xx,negmaxfn, control=list(trace=1, maxit=10, maxfeval=2000))
print(ansnegmax)

