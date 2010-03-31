# maxjn-test built on maxtestj.R
# check that optimx runs maximize for all methods

rm(list=ls())
library(optimx)

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


# source("/home/john/R-optimtest/optimx/spgfix/spg_fix.txt")


x0<-rep(pi,4)
ans.mx<-optimx(x0,maxfn,control=list(maximize=TRUE,all.methods=TRUE,save.failures=TRUE,trace=TRUE))

optansout(ans.mx, filename="./ansmx.txt")

#ans.b<-optimx(x0,maxfn,method='bobyqa',control=list(maximize=TRUE,save.failures=TRUE,trace=TRUE))
#ans.n<-optimx(x0,maxfn,method='newuoa',control=list(maximize=TRUE,save.failures=TRUE,trace=TRUE))
#ans.u<-optimx(x0,maxfn,method='uobyqa',control=list(maximize=TRUE,save.failures=TRUE,trace=TRUE))

#ans.nd<-newuoa(x0,negmaxfn,control=list(iprint=2))
#ans.ud<-uobyqa(x0,negmaxfn,control=list(iprint=2))
#ans.bd<-bobyqa(x0,negmaxfn,control=list(iprint=2))

ans.rc<-optimx(x0,maxfn,method='Rcgmin',control=list(maximize=TRUE,save.failures=TRUE,trace=TRUE))
optansout(ans.mx, filename="./ansrc.txt")
ans.rcd<-Rcgmin(x0,negmaxfn,control=list(trace=2))
optansout(ans.mx, filename="./ansrcd.txt")
# ans.rvd<-Rvmmin(x0,negmaxfn,control=list(trace=2))


x00<-c(1,2,3,4)
# Test if things work when we provide the solution!

ans.mx0<-optimx(x0,maxfn,control=list(all.methods=TRUE,maximize=TRUE,save.failures=TRUE,trace=TRUE))
optansout(ans.mx, filename="./ansmx0.txt")

#ans.spgm<-spg(x0, maxfn,control=list(maximize=TRUE,trace=TRUE))


