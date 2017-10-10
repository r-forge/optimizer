## Trig function 
# ref??

trig.f <- function(x){
  res <- trig.res(x)
  sum(res*res)
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

require("minpack.lm")
require("nlmrt")
require(optimx)
st<-c(1,2,3)
lo<-c(0.25, 0.25, 0.25)
up<-c(3,3,3)
ilo<-rep(-Inf, 3)
iup<-rep(Inf,3)
cat("Unconstrained tests\n")
anlfb1<-nlfb(st, trig.res, trig.jac, trace=TRUE)
print(anlfb1)
tmp<-readline("\ncontinue")
aopx1<-optimx(st, trig.f, trig.g, control=list(all.methods=TRUE, trace=1))
print(aopx1)
tmp<-readline("\ncontinue")
amin1<-nls.lm(st, lower=ilo, upper=iup, trig.res, trig.jac, control=nls.lm.control(nprint=1))
print(amin1)
tmp<-readline("\ncontinue")
anls1<-try(nls(st, trig.res, trig.jac, algorithm='port', trace=TRUE))
print(anls1)
tmp<-readline("\ncontinue")

cat("Bounds constrained to bounds \n")
cat("lower:")
print(lo)
cat("upper:")
print(up)



anlfb2<-nlfb(st, trig.res, lower=lo, upper=up, trig.jac, trace=TRUE)
print(anlfb2)
tmp<-readline("\ncontinue")
aopx2<-optimx(st, trig.f, trig.g, lower=lo, upper=up, control=list(all.methods=TRUE, trace=1))
print(aopx2)
tmp<-readline("\ncontinue")
amin2<-nls.lm(st, lower=lo, upper=up, trig.res, trig.jac, control=nls.lm.control(nprint=1))
print(amin2)
tmp<-readline("\ncontinue")
anls2<-try(nls(st, trig.res, trig.jac, algorithm='port', trace=TRUE, lower=lo, upper=up))
print(anls2)


lo<-rep(-10,3)
up<-rep(.4,3)
st<-lo

cat("Bounds constrained to bounds \n")
cat("lower:")
print(lo)
cat("upper:")
print(up)

anlfb3<-try(nlfb(st, trig.res, lower=lo, upper=up, trig.jac, trace=TRUE))
print(anlfb3)
tmp<-readline("\ncontinue")
aopx3<-optimx(st, trig.f, trig.g, lower=lo, upper=up, control=list(all.methods=TRUE, trace=1))
print(aopx3)
tmp<-readline("\ncontinue")
amin3<-nls.lm(st, lower=lo, upper=up, trig.res, trig.jac, control=nls.lm.control(nprint=1))
print(amin3)
tmp<-readline("\ncontinue")
anls3<-try(nls(st, trig.res, trig.jac, algorithm='port', trace=TRUE, lower=lo, upper=up))
print(anls3)

