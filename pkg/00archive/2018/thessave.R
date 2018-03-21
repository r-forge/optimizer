cat("Test Hessian save\n")

require(numDeriv)
require(optimx)

jones<-function(xx){
  x<-xx[1]
  y<-xx[2]
  ff<-sin(x*x/2 - y*y/4)*cos(2*x-exp(y))
  ff<- -ff
}

jonesg <- function(xx) {
  x<-xx[1]
  y<-xx[2]
  gx <-  cos(x * x/2 - y * y/4) * ((x + x)/2) * cos(2 * x - exp(y)) - 
    sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * 2)
  gy <- sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * exp(y)) - cos(x * 
                                                                        x/2 - y * y/4) * ((y + y)/4) * cos(2 * x - exp(y))
  gg <- - c(gx, gy)
}

xx <- c(1, 2)

ans  <- optimr(xx, jones, jonesg, method="nlminb", hessian=TRUE, control=list(trace=1))
ans
xstar <- ans$par
hn <- hessian(jones, xstar)
hn



