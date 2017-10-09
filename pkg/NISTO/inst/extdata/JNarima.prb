rm(list=ls())
#JN Mod of Ravi's arima.r test

# maximum likelihood for ARMAX model
# log likelihood - Frances & van Oest eq (3)
# http://publishing.eur.nl/ir/repub/asset/1190/ei200407.pdf

library("numDeriv")


negll <- function(inp, y, x) {
	beta <- inp[1]
	lambda <- inp[2]
        sigma <- inp[3]
	e <- 0 * y
	n <- length(e)
	for(i in 2:n) e[i] <- y[i] - beta * x[i] - lambda * y[i-1] + lambda * e[i-1]
	res<- - (n - 1)/2 * (log(2 * pi) + log(sigma^2)) - sum(e*e) / (2 *(sigma^2))
#        cat("e[n]=",e[n],"\n")
        return(-res)
}


ll <- function(inp, yy, xx) {
	beta <- inp[1]
	lambda <- inp[2]
        sigma <- inp[3]
	e <- 0 * yy
	n <- length(e)
	for(i in 2:n) e[i] <- yy[i] - beta * xx[i] - lambda * yy[i-1] + lambda * e[i-1]
	res<- - (n - 1)/2 * (log(2 * pi) + log(sigma^2)) - sum(e*e) / (2 *(sigma^2))
#        cat("e[n]=",e[n],"\n")
        return(res)
}

# grll <- function(fun, inp, yy, xx) {
grll <- function(inp, yy, xx) {
	cat("inp:")
	print(inp)
	cat("yy:")
	print(yy)
	cat("xx:")
	print(xx)
#	gr <- grad(fun, inp, yy = yy, xx = xx)
	gr <- grad(ll, inp, yy = yy, xx = xx)
#	gr <- grad(fun, inp,  method="Richardson", method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine $double.eps/7e-7), r=4, v=2, show.details=TRUE), y = y, x = x)
#	gr <- grad(fun, inp,  method="Richardson", method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine $double.eps/7e-7), r=4, v=2, show.details=TRUE), yy = yy, xx = xx)
	return(gr)
}

# generate test series, y, using predictor x and indicated beta and lambda
set.seed(123)
beta <- 0.5
lambda <- 0.95
sigma <- 100.0
n <- 50
x <- rnorm(n, 0, 10)
e <- rnorm(n, 0, sigma)
y <- 0 * e

for(i in 2:n) y[i] <- lambda * y[i-1] + beta * x[i] + e[i] - lambda * e[i-1]

library(optimx)
init <- c(1, 1, 1)

fe0<-ll(init,yy=y, xx=x)
cat("fe0: ")
print(fe0)

# gfe0<-grll(ll, init, y=y , x=x)
#gfe0<-grll(ll, init, yy=y, xx=x)
gfe0<-grll(init, yy=y, xx=x)
cat("gradient: ")
print(gfe0)


tmp<-readline("Continue?")

# fm <- optimx(init, ll, lower=c(-Inf, 0, 0), control = list(maximize=TRUE, all.methods=TRUE, trace=FALSE), yy = y, xx = x)
# Results of L-BFGS-B are displayed even though trace=FALSE

# print(fm)
# Rvmmin and Rcgmin give negative values for second parameter

# fm1 <- optimx(init, ll, lower=c(-Inf, 0, 0), control = list(maximize=TRUE, all.methods=TRUE, trace=2), yy = y, xx = x)
# print(fm1)


rvfm <- Rvmmin(init, ll, lower=c(-Inf, 0, 0), control = list(maximize=TRUE, trace=3), yy = y, xx = x)
cat("rvfm:\n")
print(rvfm)

tmp<-readline("Continue?")
rvfmg <- Rvmmin(init, ll, gr=grll, lower=c(-Inf, 0, 0), control = list(maximize=TRUE, trace=3), yy = y, xx = x)
cat("rvfmg:\n")
print(rvfmg)

tmp<-readline("Continue?")

rvfm1 <- Rvmmin(init, ll, lower=c(-1E20, 0, 0), control = list(maximize=TRUE, trace=3), yy = y, xx = x)
cat("rvfm1:\n")
print(rvfm1)

tmp<-readline("Continue?")

rvfmg1 <- Rvmmin(init, ll, gr=grll, lower=c(-Inf, 0, 0), control = list(maximize=TRUE, trace=3), yy = y, xx = x)
cat("rvfmg1:\n")
print(rvfmg1)


