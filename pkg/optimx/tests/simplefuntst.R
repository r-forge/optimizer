# simplefuntst.R
##  author: John C. Nash
rm(list=ls())
library(optimx)
sessionInfo()

# Simple Test Function 1:
tryfun.f = function(x) { 
     fun <- sum(x^2 )
#	print(c(x = x, fun = fun))
     fun
}
tryfun.g = function(x) { 
     grad<-2.0*x
     grad
}
tryfun.h = function(x) { 
     n<-length(x)
     t<-rep(2.0,n)
     hess<-diag(t)
}

strt <- c(1,2,3)
ansfgh <- optimr(strt, tryfun.f, tryfun.g, tryfun.h, method="nlm", 
     hessian=TRUE, control=list(trace=2))
ansfgh

ansall <- opm(strt, tryfun.f, tryfun.g, tryfun.h, method="ALL")
summary(ansall, order=value)
