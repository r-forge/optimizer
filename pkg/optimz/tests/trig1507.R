## Trig function 
# ref??

trig.f <- function(x){
  res <- trig.res(x)
  f <- sum(res*res)
#  cat("FV=",f," at ")
#  print(x)
  f
}

trig.res <- function(x){
   n <- length(x)
   i <- 1:n
   res <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
   return(res)
}

trig.jac <- function(x) { # not vectorized. Can it be?
## stop("Not defined")
   n <- length(x)
   J<-matrix(0,n,n)
   for (i in 1:n) {
      for (j in 1:n) {
         J[i,j]<-sin(x[j]) # we overwrite J[i,i]
      }
      J[i,i] <- (1+i) * sin(x[i])  - cos(x[i])
   }
   return(J)
}


trig.g <- function(x) { # unvectorized
  n<-length(x)
  res<-trig.res(x)
  J<-trig.jac(x)
  g<- as.vector(2.0 * ( t(J) %*% res ))
  return(g)
}


require(optest)
x<-rep(4,2)
cat("opm\n")
opt2<-optimr(x, trig.f, trig.g, method="BFGS")
opt2
opt2r<-optimr(x, trig.f, trig.g, method="Rvmmin")
opt2r
cat("====================")
x<-rep(4,4)
cat("optim(BFGS) vs Rvmmin\n")
opt4<-optim(x, trig.f, trig.g, method="BFGS")
opt4
opt4r<-optimr(x, trig.f, trig.g, method="Rvmmin")
opt4r
cat("====================")
x<-rep(4,8)
cat("optim(BFGS) vs Rvmmin\n")
opt8<-optim(x, trig.f, trig.g, method="BFGS")
opt8
opt8r<-optimr(x, trig.f, trig.g, method="Rvmmin")
opt8r

# this does NOT compute hessian in optimx -- why? 131022


# tansxu<-optimx(st, trig.f, trig.g, method="all", hessian=TRUE, control=list(kkt=TRUE, trace=1))
