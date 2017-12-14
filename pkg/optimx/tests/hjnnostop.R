library(optimr)
maxfn<-function(x) {
  n <- length(x)
  ss <- seq(1,n)
  f <- 10-(crossprod(x-ss))^2
  f <- as.numeric(f)
  return(f)
}

negmaxfn<-function(x) {
  f<-(-1)*maxfn(x)
  return(f)
}


cat("WARNING -- this example does NOT appear to terminate\n")
cat("test that maximize=TRUE works correctly\n")
# 160706 -- not set up to maximize yet, except through optimr perhaps
n<-6
xx<-rep(1,n)
ansmax<-try(hjn(xx,maxfn, control=list(maximize=TRUE,trace=1, maxit=10, maxfeval=2000)))
print(ansmax)

cat("using the negmax function should give same parameters\n")
ansnegmax<-hjn(xx,negmaxfn, control=list(trace=2))
print(ansnegmax)
