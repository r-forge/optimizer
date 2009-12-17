###################################################
## optxvmmix -- examples of use of optimx function
library(optimx)

if(!require("BB"))    stop("this requires package BB.")
if(!require(numDeriv))stop("this requires package numDeriv.")
if(!require("Rcgmin"))stop("this requires Rcgmin.")
if(!require("minqa"))stop("this requires minqa.")
if(!require("ucminf"))stop("this requires ucminf.")


###############################################################

dvm <- function (theta, mu, kappa) 
{
    1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * 
        (exp(cos(theta - mu) - 1))^kappa
}

##
rmixedvm <- function (n, mu1, mu2, kappa1, kappa2, p) {
temp <- runif(n)
n1 <- sum(temp <= p)
y <- c(rvm(n1,mu1,kappa1),rvm(n-n1,mu2,kappa2))
return(y)
}

##
rvm <- function (n, mean, k) 
{
    vm <- c(1:n)
    a <- 1 + (1 + 4 * (k^2))^0.5
    b <- (a - (2 * a)^0.5)/(2 * k)
    r <- (1 + b^2)/(2 * b)
    obs <- 1
    while (obs <= n) {
        U1 <- runif(1, 0, 1)
        z <- cos(pi * U1)
        f <- (1 + r * z)/(r + z)
        c <- k * (r - f)
        U2 <- runif(1, 0, 1)
        if (c * (2 - c) - U2 > 0) {
            U3 <- runif(1, 0, 1)
            vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
            vm[obs] <- vm[obs]%%(2 * pi)
            obs <- obs + 1
        }
        else {
            if (log(c/U2) + 1 - c >= 0) {
                U3 <- runif(1, 0, 1)
                vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
                vm[obs] <- vm[obs]%%(2 * pi)
                obs <- obs + 1
            }
        }
    }
    vm
}

#
vmmix.loglik <- function(x,y){
p <- x
value<-(-1)*sum(log(p[5]*dvm(y,p[1],p[2])+(1-p[5])*dvm(y,p[3],p[4])))
return(value)
}

y <- rmixedvm(n=500, mu1=pi/2, mu2=3*pi/2, kappa1=1.9, kappa2=2.2, p=0.67)
p <- c(pi/4,2,pi,1,0.5)

lo <- rep(0.001,5)
hi <- c(Inf, Inf, Inf, Inf, 0.999)

p0 <- c(runif(5,c(0,0.1,0,0.1,0.2),c(2*pi,5,2*pi,5,0.8)))

ans.vmmix <- optimx(par=p0, fn=vmmix.loglik, y=y, lower=lo, upper=hi, control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem vmmix\n")
print(ans.vmmix)

