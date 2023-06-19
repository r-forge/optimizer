# testbmstep.R -- see if we can run start tests

# Simple Test Function 1:
simfun.f = function(x) { 
  fun <- sum(x^2 )
  #	print(c(x = x, fun = fun))
  fun
}
simfun.g = function(x) { 
  grad<-2.0*x
  grad
}
simfun.h = function(x) { 
  n<-length(x)
  t<-rep(2.0,n)
  hess<-diag(t)
}



library(optimx)
sessionInfo()
strt<-c(1,2,3,4)
lo<-c(.1, .2, .3, .4)
up<-c(2, 3, 4, 5)
srch<-simfun.g(strt)
cat("srch:")
print(srch)

# sink(file="~/temp/tbm.txt", split=TRUE)
t1 <- try(bmstep(strt, srch, lower=lo, upper=up, bdmsk=NULL, trace=1))
cat("t1:")
print(t1)
srch<- (-1)*srch
t2 <- try(bmstep(strt, srch, lower=lo, upper=up, bdmsk=NULL, trace=1))
cat("t2:")
print(t2)
# sink()
