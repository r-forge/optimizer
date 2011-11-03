# test to run funcheck on tryfun functions
require(optfntools)

# Simple Test Function 1:
tryfun.f = function(x, scale=1) { 
     fun <- sum(x^2 )*sqrt(scale)
## if (trace) ... to be fixed
#	print(c(x = x, fun = fun))
     fun
}

tryfun.g = function(x, scale=1) { 
     grad<-2.0*x*sqrt(scale)
     grad
}

tryfun.h = function(x, scale=1) { 
     n<-length(x)
     t<-rep(2.0,n)
     hess<-diag(t)*sqrt(scale)
}




ansfc<-fnchk(xpar = rep(10.0,10), ffn=tryfun.f,scale=.33 ) 
# ?? do we want to do anything with the answer
ansfc

ansgc<-grchk(xpar=rep(10, 10), ffn=tryfun.f, ggr=tryfun.g, trace=1)
ansgc

anshc<-hesschk(xpar=rep(10,10), ffn=tryfun.f, ggr=tryfun.g, hhess=tryfun.h, trace=1)
anshc