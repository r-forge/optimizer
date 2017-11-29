require(optimr)
maxfn<-function(x) {
  n<-length(x)
  ss<-seq(1,n)
  f<-10-(crossprod(x-ss))^2
  f<-as.numeric(f)
  return(f)
}


negmaxfn<-function(x) {
  f<-(-1)*maxfn(x)
  return(f)
}

cat("test that maximize=TRUE works correctly\n")

n<-6
xx<-rep(1,n)
ansmax<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE,trace=1))
print(ansmax)

cat("using the negmax function should give same parameters\n")
ansnegmax<-optimr(xx,negmaxfn, gr="grfwd",  method="Rvmmin", control=list(trace=1))
print(ansnegmax)
