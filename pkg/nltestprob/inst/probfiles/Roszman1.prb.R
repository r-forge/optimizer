## @knitr ##Roszman1.prb
# This is file ##Roszman1.prb
rm(list=ls())
probname <- "##Roszman1"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Roszman1          (Roszman1.dat)

File Format:   ASCII
Starting Values   (lines 41 to 44)
Certified Values  (lines 41 to 49)
Data              (lines 61 to 85)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study involving
quantum defects in iodine atoms.  The response
variable is the number of quantum defects, and the
predictor variable is the excited energy state.
The argument to the ARCTAN function is in radians.

Reference:     Roszman, L., NIST (19??).  
Quantum Defects for Sulfur I Atom.






Data:          1 Response  (y = quantum defect)
1 Predictor (x = excited state energy)
25 Observations
Average Level of Difficulty
Observed Data

Model:         Miscellaneous Class
4 Parameters (b1 to b4)

pi = 3.141592653589793238462643383279E0
y =  b1 - b2*x - arctan[b3/(x-b4)]/pi  +  e


Starting Values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =      0.1         0.2         2.0196866396E-01  1.9172666023E-02
b2 =     -0.00001    -0.000005   -6.1953516256E-06  3.2058931691E-06
b3 =   1000        1200           1.2044556708E+03  7.4050983057E+01
b4 =   -100        -150          -1.8134269537E+02  4.9573513849E+01

Residual Sum of Squares:                    4.9484847331E-04
Residual Standard Deviation:                4.8542984060E-03
Degrees of Freedom:                                 21
Number of Observations:                             25










Data:   y           x
0.252429    -4868.68
0.252141    -4868.09
0.251809    -4867.41
0.297989    -3375.19
0.296257    -3373.14
0.295319    -3372.03
0.339603    -2473.74
0.337731    -2472.35
0.333820    -2469.45
0.389510    -1894.65
0.386998    -1893.40
0.438864    -1497.24
0.434887    -1495.85
0.427893    -1493.41
0.471568    -1208.68
0.461699    -1206.18
0.461144    -1206.04
0.513532     -997.92
0.506641     -996.61
0.505062     -996.31
0.535648     -834.94
0.533726     -834.66
0.568064     -710.03
0.612886     -530.16
0.624169     -464.17

"

#- Note: environment / list "counters" must already exist

if (exists("pe")) { 
  rm("pe")  
}

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

#- nls format expression
Roszman1.formula <- ( y ~ b1 - b2*x - atan(b3/(x-b4))/pi )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Roszman1")))
# Optimization test function roszman1
# roszman1 from NISTnls
# ??ref...


roszman1.f <- function(x) {
   res<-roszman1.res(x)
   f<-sum(res*res)
}

roszman1.res <- function(b) {
   xx<-Roszman1$x
   yy<-Roszman1$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<- b1 - b2*xx - atan(b3/(xx-b4))/pi  - yy
   return(res)
}

# roszman1 - Jacobian
roszman1.jac <- function(b) {

   xx<-Roszman1$x
   yy<-Roszman1$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   J<-matrix(0,m,n) # define the size of the Jacobian
   J[,1] <- 1
   J[,2] <- -xx
   tmp <- 1/(pi*(1+(b3/(xx-b4))^2))
   J[,3] <- tmp /(xx - b4)
   J[,4] <-  tmp * b3 / (xx - b4)^2 
   
   J
}

roszman1.h <- function(x) {
stop("not defined")
   JJ<-roszman1.jac(x)
   H <- t(JJ) %*% JJ
   res<-roszman1.res(x)
stop("not defined")

}

roszman1.g<-function(x) {
#   stop("not defined")
   JJ<-roszman1.jac(x)
   res<-roszman1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

roszman1.fgh<-function(x) {
   f<-roszman1.f(x)
   g<-roszman1.g(x)
   H<-roszman1.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}


start0<-rep(1,4)
names(start0) <- c("b1","b2","b3","b4")
start1<-c(0.1, -0.00001, 1000, -100)
names(start1) <- c("b1","b2","b3","b4")
start2<-c(0.2, -0.000005, 1200, -150)
names(start2) <- c("b1","b2","b3","b4")



NISTRoszman1<-list()
NISTRoszman1$value<- 4.9484847331E-04
NISTRoszman1$par<-c( 2.0196866396E-01,-6.1953516256E-06, 1.2044556708E+03, -1.8134269537E+02)
NISTRoszman1$ses<-c( 1.9172666023E-02, 3.2058931691E-06, 7.4050983057E+01, 4.9573513849E+01)

library(numDeriv)

cat("f(start0) =", roszman1.f(start0),"\n")
JA <- roszman1.jac(start0)
JN <- jacobian(roszman1.res, start0)
cat("max(abs(diff)) =",max(abs(JA - JN)),"\n")


   
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Roszman1)
Try(fm1 <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
               start = c(b1 = 0.1, b2 = -0.00001, b3 = 1000, b4 = -100),
               trace = TRUE))
Try(fm1a <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
                start = c(b1 = 0.1, b2 = -0.00001, b3 = 1000, b4 = -100),
                alg = "port", trace = TRUE))
Try(fm2 <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
               start = c(b1 = 0.2, b2 = -0.0000015, b3 = 1200, b4 = -150),
               trace = TRUE))
Try(fm2a <- nls(y ~ b1 - b2*x - atan(b3/(x-b4))/pi, data = Roszman1,
                start = c(b1 = 0.2, b2 = -0.0000015, b3 = 1200, b4 = -150),
                alg = "port", trace = TRUE))

library(nlsr)

Roszman1nlsr0 <- nlxb(start=start0, formula=Roszman1.formula, data=mypdata, trace=TRUE)
print(Roszman1nlsr0)

Roszman1nlsr1 <- nlxb(start=start1, formula=Roszman1.formula, data=mypdata, trace=TRUE)
print(Roszman1nlsr1)

Roszman1nlsr2 <- nlxb(start=start2, formula=Roszman1.formula, data=mypdata, trace=TRUE)
print(Roszman1nlsr2)

library(optimr)

Roszman1opm0fwd <- opm(start0, roszman1.f, "grfwd", method="ALL")
summary(Roszman1opm0fwd, order=value)

Roszman1opm1fwd <- opm(start1, roszman1.f, "grfwd", method="ALL")
summary(Roszman1opm1fwd, order=value)

Roszman1opm2fwd <- opm(start2, roszman1.f, "grfwd", method="ALL")
summary(Roszman1opm2fwd, order=value)

Roszman1opm0 <- opm(start0, roszman1.f, roszman1.g, method="ALL")
summary(Roszman1opm0, order=value)

Roszman1opm1 <- opm(start1, roszman1.f, roszman1.g, method="ALL")
summary(Roszman1opm1fwd, order=value)

Roszman1opm2 <- opm(start2, roszman1.f, roszman1.g, method="ALL")
summary(Roszman1opm2, order=value)

