fn <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  (x-1)^2 + (x-y-3)^4 + 0.0001*(y-4)^4 + 1
}
start <- c(x=0, y=10)
start
set.seed(12345) # REALLY IMPORTANT!
a1s <- optim(start, fn, method="SANN", control=list(maxit=20))
a1s
a1n <- optim(start, fn, method="Nelder-Mead", control=list(maxit=20))
a1n
a1sn <- optim(a1s$par, fn, method="Nelder-Mead", control=list(maxit=20))
a1sn