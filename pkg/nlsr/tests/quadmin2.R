# Test of nlsr hessian and gradient

fexp <- ~ (x-1)^2 + (y-3)^2
fexp
library(nlsr)
myder <- fnDeriv(fexp, c("x","y"), hessian=TRUE)
myder
strt <- c(x=0, y=-0)
astrt <- myder(strt)
myd <- function(prm=strt){
  x <- prm[1]
  y <- prm[2]
  val <- myder(x, y)
}
astrt <- myd(strt)
astrt
hm1 <- as.matrix(attr(astrt, "hessian"))
hm1
hm <- matrix(hm1, nrow=2)
hm
g1 <- attr(astrt, "gradient")
g1 <- as.vector(g1)
g1
delta <- solve(hm, g1)
delta
new <- strt-delta
new
newf<-myd(new)
newf
