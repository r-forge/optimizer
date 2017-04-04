## Trig function 
# ref  More' Garbow and Hillstrom, 1981, problem 26, from
#  Spedicato (their ref. 25)

# trig.f <- function(x){
#  res <- trig.res(x)
#  f <- sum(res*res)
#  cat("FV=",f," at ")
#  print(x)
#  f
# }

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
   attr(J, "gradient") <- J
   J
}


# trig.g <- function(x) { # unvectorized
#  n<-length(x)
#  res<-trig.res(x)
#  J<-trig.jac(x)
#  g<- as.vector(2.0 * ( t(J) %*% res ))
#  return(g)
#}
