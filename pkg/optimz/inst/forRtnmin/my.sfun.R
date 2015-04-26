my.sfun <- function(x){
##------------------------------------------
## Get function value and gradient for the
## objective function.
##------------------------------------------
## Usage: [f,g] = my_sfun(x)
##------------------------------------------

## Specify matrix and right-hand side for quadratic function

   n <- length(x) 
   A <- diag(1:n) 
   for (i in 2:n) {
     A[i,i-1] <- 1
     A[i-1,i] <- 1
   }
   b <- rep(1, n)

## Evaluate function and gradient
   g <- as.numeric(crossprod(A, x) - b)
   f  <- as.numeric(0.5*crossprod(x, g-b))
##   print(A)
   list(f=f, g=g)
}

fonly.sfun<-function(x){
   my.sfun(x)$f
}

gonly.sfun<-function(x){
   my.sfun(x)$g
}