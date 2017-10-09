## @knitr ##Misra1a.prb
# This is file ##Misra1a.prb
rm(list=ls())
probname <- "##Misra1a"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Misra1a           (Misra1a.dat)

File Format:   ASCII
Starting Values   (lines 41 to 42)
Certified Values  (lines 41 to 47)
Data              (lines 61 to 74)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study regarding
dental research in monomolecular adsorption.  The
response variable is volume, and the predictor
variable is pressure.

Reference:     Misra, D., NIST (1978).  
Dental Research Monomolecular Adsorption Study.







Data:          1 Response Variable  (y = volume)
1 Predictor Variable (x = pressure)
14 Observations
Lower Level of Difficulty
Observed Data

Model:         Exponential Class
2 Parameters (b1 and b2)

y = b1*(1-exp[-b2*x])  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   500         250           2.3894212918E+02  2.7070075241E+00
b2 =     0.0001      0.0005      5.5015643181E-04  7.2668688436E-06

Residual Sum of Squares:                    1.2455138894E-01
Residual Standard Deviation:                1.0187876330E-01
Degrees of Freedom:                                12
Number of Observations:                            14












Data:   y               x
10.07E0      77.6E0
14.73E0     114.9E0
17.94E0     141.1E0
23.93E0     190.8E0
29.61E0     239.9E0
35.18E0     289.0E0
40.02E0     332.8E0
44.82E0     378.4E0
50.76E0     434.8E0
55.05E0     477.3E0
61.01E0     536.8E0
66.40E0     593.1E0
75.47E0     689.1E0
81.78E0     760.0E0

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
Misra1a.formula <- ( y ~ b1*(1-exp(-b2*x)) )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Misra1a")))
# Optimization test function misra1a
# misra1a from NISTnls
# ??ref...


misra1a.f <- function(x) {
   res<-misra1a.res(x)
   f<-sum(res*res)
}

misra1a.res <- function(b) {
   xx<-Misra1a$x
   yy<-Misra1a$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   res<-b1*(1-exp(-b2*xx)) - yy
   return(res)
}

# misra1a - Jacobian
misra1a.jac <- function(b) {
  xx<-Misra1a$x
  yy<-Misra1a$y
  n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   J[,1] <- (1-exp(-b2*xx))
   J[,2] <- b1*xx*exp(-b2*xx)
   J 
}

misra1a.h <- function(x) {
stop("not defined")
   JJ<-misra1a.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1a.res(x)
stop("not defined")

}

misra1a.g<-function(x) {
#   stop("not defined")
   JJ<-misra1a.jac(x)
   res<-misra1a.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1a.fgh<-function(x) {
   f<-misra1a.f(x)
   g<-misra1a.g(x)
   H<-misra1a.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1<-c(500, 0.0001)
names(start1) <- c("b1", "b2")

start2<-c(250, 0.0005)
names(start2) <- c("b1", "b2")

## Examples

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Misra1a)
Try(fm1 <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
               start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
                alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
               start = c(b1 = 250, b2 = 0.0005) ))
Try(fm2a <- nls(y ~ b1*(1-exp(-b2*x)), data = Misra1a, trace = TRUE,
                alg = "port", start = c(b1 = 250, b2 = 0.0005) ))
Try(fm3 <- nls(y ~ 1-exp(-b2*x), data = Misra1a, trace = TRUE,
               start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ 1-exp(-b2*x), data = Misra1a, trace = TRUE,
               start = c(b2 = 0.0005), algorithm = "plinear" ))

## Using a self-starting model
Try(fm5 <- nls(y ~ SSasympOrig(x, Asym, lrc), data = Misra1a))

# nlsr
library(nlsr)

misra1anlxb1 <- nlxb(start=start1, formula=Misra1a.formula, data=mypdata, trace=TRUE)
print(misra1anlxb1)


misra1anls2 <- nls(start=start2, formula=Misra1a.formula, data=mypdata, trace=TRUE)
print(misra1anls2)

misra1anlxb2 <- nlxb(start=start2, formula=Misra1a.formula, data=mypdata, trace=TRUE)
print(misra1anlxb2)


# optimr
library(optimr)