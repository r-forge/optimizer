 four<-function(xt) {
     n<-length(xt)
     res<-xt-rep(4,n)
     f<- exp(-sum(res*res)/10)
     f
 }

