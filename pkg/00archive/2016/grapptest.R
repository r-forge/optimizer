# Test for gradient approx in optimr and opm
nondia.f <- function(x){
   n <- length(x)
   if (n < 2) { stop("n too small for NONDIA function") }
   x1 <- x[1]
   y <- x^2 - x1
   z <- x - 1
   val <- sum(4 * y^2 + z^2)
}

nondia.g <- function(x){
   n <- length(x)
   if (n < 2) { stop("n too small for NONDIA function") }
   x1 <- x[1]
   y <- x^2 - x1
   z <- x - 1
   g <- rep(NA,n)
   g1 <- 16 * y * x
   g1[1] <- g1[1] - 8 * sum(y)
   g2 <- 2 * z
   g <- g1 + g2
}

tmp <- readline("continue")

library(optimrx)

s4 <- rep(-1,4)

sol4 <- opm(s4, nondia.f, nondia.g, method="ALL")
summary(sol4, order=value)
tmp <- readline("continue")

sol4f <-opm(s4, nondia.f, gr="grfwd", method="ALL")
summary(sol4f, order=value)
tmp <- readline("continue")
 
sol4c <-opm(s4, nondia.f, gr="grcentral", method="ALL")
summary(sol4c, order=value)
tmp <- readline("continue")

sol4b <-opm(s4, nondia.f, gr="grback", method="ALL")
summary(sol4b, order=value)
tmp <- readline("continue")

sol4n <-opm(s4, nondia.f, gr="grnd", method="ALL")
summary(sol4n, order=value)
tmp <- readline("continue")

sol4dn <-opm(s4, nondia.f, gr=NULL, method="ALL")
summary(sol4dn, order=value)
tmp <- readline("continue")

bfgs4n <- optimr(s4, nondia.f, gr="grnd", method="BFGS")
print(bfgs4n)

Rcg4n <- optimr(s4, nondia.f, gr="grnd", method="Rcgmin")
print(Rcg4n)

Rcg4 <- optimr(s4, nondia.f, gr=nondia.g, method="Rcgmin")
print(Rcg4)


Rvm4n <- optimr(s4, nondia.f, gr="grnd", method="Rvmmin")
print(Rvm4n)

Rvm4 <- optimr(s4, nondia.f, gr=nondia.g, method="Rvmmin")
print(Rvm4)
tmp <- readline("continue")
 

sol4fu <-opmu(s4, nondia.f, gr="grfwd", method="ALL")
summary(sol4fu, order=value)
tmp <- readline("continue")
 
sol4nu <-opmu(s4, nondia.f, gr="grnd", method="ALL")
summary(sol4nu, order=value)
tmp <- readline("continue")
 
