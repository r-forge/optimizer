rm(list=ls())
# fname<-paste(format(Sys.time(), "%Y%m%d%H%M"),"-btRvmmin.out",sep='')
# sink(fname, append=TRUE, split=TRUE)
require("optimz")
#####################

# This test script illustrates the use of bounds in optimr() with the
# optimizers Rvmmin and L-BFGS-B, as well as a Kuhn Karush Tucker check 
# on the final parameters from the second optimization.


# Simple bounds test for n=4
bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

n<-4
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
# bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
    lower[i]<-1.0*(i-1)*(n-1)/n
    upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)

cat("Rvmmin \n\n")

abtrvm <- optimr(xx, bt.f, bt.g, lower, upper, method="Rvmmin", control=list(trace=4))
print(abtrvm)

alb<-optimr(xx,bt.f, bt.g, lower=lower, upper=upper, method="L-BFGS-B", control=list(trace=3))
print(alb)

alkkt <- kktchk(alb$par, bt.f, bt.g, hess=NULL, upper=upper, lower=lower,  maxfn=FALSE, control=list(trace=1))
print(alkkt)


#sink()

