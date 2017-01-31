##
##  Test functions with many local minima
##


#-- 2-dim. test functions --------------
# No. of Vars.:  n = 2
# Bounds:  -1 <= xi <= 1, i=1, 2
# Local minima:  
# Minimum:  
# Solution:  
fnTrefethen <- function(x) {
    x1 <- x[1]; x2 <- x[2]
    return(exp(sin(50*x1)) + sin(60*exp(x2)) + sin(70*sin(x1)) +
               sin(sin(80*x2)) - sin(10*(x1+x2)) + (x1^2 + x2^2)/4)
}


#-- 3-dim. test functions --------------
# No. of Vars.:  n = 3
# Bounds:  -1 <= xi <= 1, i=1..3
# Local minima:  
# Minimum:  
# Solution:  
fnWagon <- function(x) {
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]
    return(exp(sin(50*x1)) + sin(60*exp(x2))*sin(60*x3) +
               sin(70*sin(x1))*cos(10*x3) +
               sin(sin(80*x2)) - sin(10*(x1+x3)) + (x1^2 + x2^2 + x3^2)/4)
}


#-- 4-dim. test functions --------------
##  Colville function
# No. of Vars.:  n = 4
# Bounds:  -10 <= xi <= 10, i=1..4
# Local minima:  
# Minimum:  0
# Solution:  (1.0, 1.0, 1.0, 1.0)
fnColville <- function(x) {
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
    # return
    100.0*(x1^2 - x2)^2 + (x1 - 1)^2 + (x3 - 1)^2 + 90*(x3^2 - x4)^2 +
        10.1*((x2 - 1)^2 + (x4 - 1)^2) + 19.8*(x2 - 1)*(x4 - 1)
}


#-- 5-dim. test functions --------------
##  Hald-Madsen function (non-smooth minimax function !)
# No. of Vars.:  n = 5
# Bounds:  -1 <= xi <= 1, i=1..5
# Local minima:  
# Minimum:  
# Solution:  
fnHald = function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    if (all(x == 1)) return(exp(1))
    t <- -1 + (c(1:21) - 1)/10
    v <- (x[1]+x[2]*t) / (1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
    max(abs(v))
}

##   Ackley's function
# No. of Vars.:  n >= 2
# Bounds:  -30 <= xi <= 30
# Local minima:  many, 2 global minima
# Minimum:  -13.37957501 (n=5)
# Solution:  (+-1.5157285, -1.1151432, -1.1096511, -1.1038473, -0.7471183)
fnAckley <- function(x) {
    n <- length(x)
    x1 <- x[-n]; xn <- x[-1]
    sum( exp(-0.2)*sqrt(x1^2 + xn^2) + 3*(cos(2*x1) + sin(2*xn)))
}

##  Modified Langerman function
# No. of Vars.:  1 <= n <= 10
# Bounds:  0 <= xi <= 10
# Local minima:  many, irregular
# Minimum:  0.96499992 (n=5)
# Solution:  (8.074000, 8.777001, 3.467004, 1.863013, 6.707995)
fnLangerman <- function(x) {
    n <- length(x)
    a <- c(0.806, 0.517, 0.100, 0.908, 0.965)
    A <- matrix(c(
        9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020,
        9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374,
        8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982,
        2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426,
        8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567),
        nrow=5, byrow=TRUE)
    fun <- 0.0
    for (i in 1:5) {
        b <- sum((x-A[i,1:n])^2)
        fun <- fun - sum( a[i] * exp(-1/pi*b) * cos(pi*b) )
    }
    return(fun)
}


#-- 6-dim. test functions --------------
##  Hartmann's function
# No. of Vars:  n = 6
# Bounds: all xi in [0, 1]
# Local minima: 6
# Minimum: -3.322368
# Solution: (0.20169,0.150011,0.476874,0.275332,0.311652,0.65730)
fnHartmann6 <- function(x) {
    n <- length(x)
    a <- c(1.0, 1.2, 3.0, 3.2)
    A <- matrix(c(10.0,  0.05, 3.0, 17.0,
                  3.0, 10.0,  3.5,  8.0,
                  17.0, 17.0,  1.7,  0.05,
                  3.5,  0.1, 10.0, 10.0,
                  1.7,  8.0, 17.0,  0.1,
                  8.0, 14.0,  8.0, 14.0), nrow=4, ncol=6)
    B  <- matrix(c(.1312,.2329,.2348,.4047,
                   .1696,.4135,.1451,.8828,
                   .5569,.8307,.3522,.8732,
                   .0124,.3736,.2883,.5743,
                   .8283,.1004,.3047,.1091,
                   .5886,.9991,.6650,.0381), nrow=4, ncol=6)
    fun <- 0.0
    for (i in 1:4) {
        fun <- fun - a[i] * exp(-sum(A[i,]*(x-B[i,])^2))
    }
    return(fun)
}


#-- 10-dim. test functions -------------
##  Michalewicz function
# No. of Vars.:  n >= 1
# Bounds:  0 <= xi <= pi
# Local minima:  n!
# Minimum: -4.6877 (n=5); -9.66015176 (n=10)
# Solution: (2.20290555, 1.57079635, 1.28499159, 1.92305848, 1.72046978,
# (n=10)     1.57079633, 1.45441397, 1.75608651, 1.65571740, 1.57079631)
fnMichalewicz <- function(x) {
    n <- length(x); m <- 10
    -sum( sin(x) * (sin((1:n)*x^2/pi))^(2*m) )
}


#-- n-dim. test functions --------------
##  Rastrigin function
# No. of Vars.:  n >= 1
# Bounds:  -5.12 <= xi <= 5.12
# Local minima:  many
# Minimum:  0.0
# Solution:  xi = 0, i=1:n
fnRastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10*cos(2*pi*x))
}

##  Whitley's function
# No. of Vars.:  n >= 2
# Bounds:  -10 <= xi <= 10
# Local minima:  many
# Minimum:  0
# Solution:  xi = 1, i=1,n
fnWhitley <- function(x) {
    n <- length(x)
    fun <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            fun <- fun + (100*(x[i]^2-x[j])^2 + (1-x[j])^2)^2/4000 - 
                cos(100*(x[i]^2-x[j])^2 + (1-x[j])^2) +1
        }
    }
    return(fun)
}

##  ARWhead
# No. of Vars:  n >= 2
# Bounds:
# Local minima:
# Minimum:
# Solution:
fnARWhead <- function(x) {
    #  sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2;
    n <- length(x)
    xn <- x[-n]
    sum((xn^2 + x[n]^2)^2 - 4*xn + 3.0)
}

##  Rana function
# No. of Vars.:  n >= 2
# Bounds:  -520 <= xi <= 520
# Local minima:  many
# Minimum:  -512.75316243
# Solution:  xi = -514.041683, i=1..n
fnRana <- function(x) {
    n <- length(x)
    xi <- x; xj <- x[c(2:n,1)]
    return( sum(xi*sin(sqrt(abs(xj+1-xi)))*cos(sqrt(abs(xj+1+xi)))
                + (xj+1)*cos(sqrt(abs(xj+1-xi)))*sin(sqrt(abs(xj+1+xi))))/n )
}
