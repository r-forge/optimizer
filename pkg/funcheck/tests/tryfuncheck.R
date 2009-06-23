# test to run funcheck on tryfun functions
require(funcheck)
ansfc<-funcheck(xpar = rep(10.0,10), fname="tryfun") 
# ?? do we want to do anything with the answer

source("tryfun.R") # load in the functions
ansft<-funtest(xpar = rep(10.0,10), fn=tryfun.f, gr=tryfun.g, hess=tryfun.h)


