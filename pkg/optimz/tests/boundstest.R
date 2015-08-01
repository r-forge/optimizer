# Bounds Test bt.R
# ref BT.RES in Nash and Walker-Smith (1987)

bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

bt.badsetup<-function(n){
   x<-rep(0,n)
   lower<-rep(0,n)
   upper<-lower # to get arrays set
   bdmsk<-rep(1,n)
   bdmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      x[i]<-2.2*i-n
      lower[i]<-1.0*(i-1)*(n-1)/n
      upper[i]<-1.0*i*(n+1)/n
   }
   result<-list(x=x, lower=lower, upper=upper, bdmsk=bdmsk)
}

bt.setup0<-function(n){
   x<-rep(0,n)
   lower<-rep(0,n)
   upper<-lower # to get arrays set
   bdmsk<-rep(1,n)
   bdmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      lower[i]<-1.0*(i-1)*(n-1)/n
      upper[i]<-1.0*i*(n+1)/n
   }
   x<-0.5*(lower+upper)
   result<-list(x=x, lower=lower, upper=upper, bdmsk=bdmsk)
}

bt.res<-function(x){
    stop(" RESIDUAL NOT YET DEFINED ")
}


bt.jac<-function(x){
    stop(" JACOBIAN NOT YET DEFINED ")
      
}

require(optimz)

n <- 10
setok <- bt.setup0(n)
start <- setok$x
lo <- setok$lower
up <- setok$upper 
abt <- opm(start, bt.f, bt.g, lower=lo, upper=up, method="ALL") # Must be capitalized
print(summary(abt, order=value))


abtn <- opm(start, bt.f, gr="grnd", lower=lo, upper=up, method="ALL")
print(summary(abtn, order=value))

