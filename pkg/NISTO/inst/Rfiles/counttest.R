## @knitr counttest.R

rm(list=ls())
sq<-function(x){
   nn<-length(x)
   yy<-1:nn
   f<-sum((yy-x)^2)
#   cat("Fv=",f," at ")
#   print(x)
   counters$kf <- counters$kf+1
   f
}
sq.g <- function(x){
   nn<-length(x)
   yy<-1:nn
   counters$kg <- counters$kg+1
   gg<- 2*(x - yy)
}


require(optimrx)
x<-rep(pi,3)

counters <- new.env()
counters$kf <- 0
counters$kg <- 0

soln1 <- optimr(x, sq, method="Nelder-Mead")
print(soln1)
cat("counters$kf =",counters$kf,"  counters$kg=",counters$kg,"\n")
counters$kf <- 0
counters$kg <- 0
soln2 <- optimr(x, sq, sq.g, method="BFGS")
print(soln2)
cat("counters$kf =",counters$kf,"  counters$kg=",counters$kg,"\n")


