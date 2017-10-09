# test function from book by Owen Jones et al 

# function is 2 variables to be maximized
#  f(x,y)=sin(x^2/2 - y^2/4)*cos(2x-exp(y))

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


library(optimx)

xx<-0.5*c(pi,pi)
ans<-optimx(xx, jones, jonesg, method="all", control=list(trace=1))
summary(ans, order=value)

lo<-c(0,0)
up<-c(1.2*pi, 1.2*pi)
ansb<-optimx(xx, jones, jonesg, lower=lo, upper=up, control=list(all.methods=TRUE, trace=1))
summary(ansb, order=value)

ansbn<-optimx(xx, jones, lower=lo, upper=up, control=list(all.methods=TRUE, trace=1))
summary(ansbn, order=value)

