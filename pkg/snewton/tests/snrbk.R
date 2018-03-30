require(snewton)
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

x0 <- c(-1.2, 1)

xx <- x0

# sink("mbrn1-170408.txt", split=TRUE)
t1 <- snewton(x0, fn=f, gr=gr, hess=h, control=list(trace=1))
print(t1)

t1m <- snewtonm(x0, fn=f, gr=gr, hess=h, control=list(trace=1))
print(t1m)

# we can also use nlm and nlminb
fght <- function(x){
  ## combine f, g and h into single function for nlm
     ff <- f(x)
     gg <- gr(x)
     hh <- h(x)
     attr(ff, "gradient") <- gg
     attr(ff, "hessian") <- hh
     ff
}

## Seems not to work as a Newton method
t1nlm <- nlm(fght, x0, hessian=TRUE, print.level=1)
print(t1nlm)

## BUT ... it looks like nlminb is NOT using a true Newton-type method
t1nlminb <- nlminb(x0, f, gradient=gr, hessian=h, control=list(trace=1))
print(t1nlminb)
# and call them from optimx (i.e., test this gives same results)

library(optimx)
t1nlmo <- optimr(x0, f, gr, hess=h, method="nlm", control=list(trace=1))
print(t1nlmo)

## FOLLOWING SHOWS UP ERRORS??
t1nlminbo <- optimr(x0, f, gr, hess=h, method="nlminb", control=list(trace=1))
print(t1nlminb)


# sink()
