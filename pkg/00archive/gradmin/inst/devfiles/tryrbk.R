# source(file="SNewton.R")
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
cat("Initial fn=", f(xx),"\n")
xn <- rep(Inf, 2)
while(TRUE) {
   gg <- gr(xx)
   HH <- h(xx)
   st <- solve(HH, -gg)
   xn <- xx + st
   cat("New fn = ", f(xn),"\n")
   tmp <- readline("continue?")
   xx <- xn
}

# sink("mbrn1-170408.txt", split=TRUE)
sol1 <- snewton(x0, fn=f, gr=gr, hess=h, control=list(trace=2))
print(sol1)
# sink()

