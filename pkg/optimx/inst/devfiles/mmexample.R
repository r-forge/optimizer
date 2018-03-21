##' Wood function (4 arguments 'x1' ... 'x4')
fwood <- function(x1,x2,x3,x4) {
  100*(x1^2-x2)^2 + (1-x1)^2 + 90*(x3^2-x4)^2 + (1-x3)^2 +
    10.1*((1-x2)^2 + (1-x4)^2) + 19.8*(1-x2)*(1-x4)
}
## automatically construct correct gradient and hessian:
woodf.gh <- function(x) {
  stopifnot(is.numeric(x))
  woodGH <- deriv3(body(fwood)[[2]],
                   c("x1","x2","x3","x4"), function.arg=TRUE)
  if(length(x) == 4)
    woodGH(x[1],x[2],x[3],x[4])
  else if(is.matrix(x) && ncol(x) == 4)
    woodGH(x[,1], x[,2], x[,3], x[,4])
  else stop("'x' must have length 4 or be a matrix with 4 columns")
}

x1 <- c(0, 0, 0, 0)
x2 <- c(1, 1, 1, 1)
x3 <- c(1, 2, 3, 4)

woodf.gh(rbind(0, 1, 1:4))
## [1]   42.0    0.0 2514.4
## attr(,"gradient")
## x1    x2   x3     x4
## [1,]   -2 -40.0   -2  -40.0
## [2,]    0   0.0    0    0.0
## [3,] -400 279.6 5404 -819.6
## attr(,"hessian")
## , , x1

##x1   x2 x3 x4
##[1,]   2    0  0  0
##[2,] 802 -400  0  0
##[3,] 402 -400  0  0

##, , x2

##x1    x2 x3   x4
##[1,]    0 220.2  0 19.8
##[2,] -400 220.2  0 19.8
##[3,] -400 220.2  0 19.8

##, , x3

##x1 x2   x3    x4
##[1,]  0  0    2     0
##[2,]  0  0  722  -360
##[3,]  0  0 8282 -1080

##, , x4

##x1   x2    x3    x4
##[1,]  0 19.8     0 200.2
##[2,]  0 19.8  -360 200.2
##[3,]  0 19.8 -1080 200.2
