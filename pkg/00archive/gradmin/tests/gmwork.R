#Rosenbrock banana valley function
f <- function(x){
return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
#gradient
gr <- function(x){
return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
#Hessian
h <- function(x) {
a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
return(matrix(c(a11, a21, a21, 200), 2, 2))
}

fg <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  val
}
fgh <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  attr(val,"hessian") <- h(x)
  val
}

source("/home/john/rsvnall/optimizer/pkg/gradmin/R/terminate.R")
source("/home/john/rsvnall/optimizer/pkg/gradmin/R/gradminu.R")
source("/home/john/rsvnall/optimizer/pkg/gradmin/R/lsback.R")
source("/home/john/rsvnall/optimizer/pkg/gradmin/R/vmrf.R")


x0 <- c(-1.2, 1)

sdefault <- gradminu(x0, fn=f, gr=gr, hess=h, control=list(srchdirn=snewt, trace=1, watch=TRUE))
print(sdefault)

stop()
tmp <- readline("Now for the Wood problem")

#Example 2: Wood function
#
#gradient:
wood.g <- function(x){
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}
wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  return(res)
}
#hessian:
wood.h <- function(x){
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}
#################################################
x0w <- c(-3, -1, -3, -1)
wdefault <- gradminu(x0w, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=1, watch=TRUE))
print(wdefault)

