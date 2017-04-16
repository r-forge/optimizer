
# Try testing calls to see what is transferred (eventually test also ...)
# setup
x0<-c(1,2,3,4)
fnt <- function(x, fscale=10){
  yy <- length(x):1
  val <- sum((yy*x)^2)*fscale
}
grt <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  #    gg <- rep(NA,nn)
  gg <- 2*(yy^2)*x*fscale
  gg
}

hesst <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  hh <- diag(2*yy^2*fscale)
  hh
}

t1 <- gradminu(x0, fnt, grt, hesst, control=list(trace=2), fscale=3.0)
t1
