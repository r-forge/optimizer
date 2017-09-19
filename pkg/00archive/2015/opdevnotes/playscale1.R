
# playscale1.R
library(numDeriv)
f <- function(x, t) { x[1]*exp(-x[2]*t) }
t <- 5
x <- c(5, .5)
f(x,t)
s <- function(q, t, ps) {
tpar <- q * ps
val <- f(tpar, t)
}
ps <- c(1, .1)
q <- c(5,5)
s(q, t, ps)
print(s(q, t, ps))
gn <- grad(f, x, t=t)
print(gn)
ga <- function(x, t){
g1 <- exp(-x[2]*t)
g2 <- -x[1]*t*g1
gr <- c(g1, g2)
}
print(ga(x,t))
gs <- function(q, t, ps){
tpar <- q * ps
ga(tpar, t)*ps
}
print(gs(q, t, ps))
gns <- grad(s, q, ps=ps, t=t)
print(gns)
hn <- hessian(f, x, t=t)
print(hn)
hns <- hessian(s, q, ps=ps, t=t)
print(hns)
Z <- diag(ps)
Z
Z1<-diag(1/ps)
gns
gn %*% Z
cat(" Z hn Z - hns:\n")
print(Z %*% hn %*% Z - hns)
