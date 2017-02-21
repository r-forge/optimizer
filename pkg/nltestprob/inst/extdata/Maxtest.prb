# maxtest function -- a function to maximize
# ref??


maxtest.f <- function(x) {
 n <- length(x)
 res<-(seq(3,(n+2))+0.25*rep(1,n)-x)
 ff<-exp(-sum(res*res)/25)
 ff
 }

maxtest.g <- function(x) {
 n <- length(x)
 res<-seq(3,(n+2))+0.25*rep(1,n)-x
 g<- 2*res*exp(-sum(res*res)/25)/25
 g
}

n<-10
xx<-rep(0,n)
ansmax<-Rvmmin(xx,maxtest.f, maxtest.g,control=list(maximize=TRUE,trace=1))
ansmaxn<-Rvmmin(xx,maxtest.f, control=list(maximize=TRUE,trace=1))

