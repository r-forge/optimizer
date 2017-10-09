## @knitr ##Misra1c.prb
# This is file ##Misra1c.prb
rm(list=ls())
probname <- "##Misra1c"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Misra1c           (Misra1c.dat)

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
Dental Research Monomolecular Adsorption.







Data:          1 Response  (y = volume)
1 Predictor (x = pressure)
14 Observations
Average Level of Difficulty
Observed Data

Model:         Miscellaneous Class
2 Parameters (b1 and b2)

y = b1 * (1-(1+2*b2*x)**(-.5))  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   500         600           6.3642725809E+02  4.6638326572E+00
b2 =     0.0001      0.0002      2.0813627256E-04  1.7728423155E-06

Residual Sum of Squares:                    4.0966836971E-02
Residual Standard Deviation:                5.8428615257E-02
Degrees of Freedom:                                12
Number of Observations:                            14












Data:   y            x 
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
Misra1c.formula <- ( y ~ b1 * (1-(1+2*b2*x)**(-.5)) )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Misra1c")))
# Optimization test function misra1c
# misra1c from NISTnls
# ??ref...


misra1c.f <- function(x) {
   res<-misra1c.res(x)
   f<-sum(res*res)
}

misra1c.res <- function(b) {
   xx<-Misra1c$x
   yy<-Misra1c$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1 * (1-(1+2*b2*xx)**(-.5)) - yy
   return(res)
}

# misra1c - Jacobian

misra1c.jac <- function(b) {
# stop("not defined")
   xx<-Misra1c$x
   yy<-Misra1c$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   J[,1] <- (1-(1+2*b2*xx)**(-.5))
   J[,2] <- b1 * xx * ((1+2*b2*xx)**(-1.5)) 
   J
}

misra1c.h <- function(x) {
stop("not defined")
   JJ<-misra1c.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1c.res(x)
stop("not defined")

}

misra1c.g<-function(x) {
#   stop("not defined")
   JJ<-misra1c.jac(x)
   res<-misra1c.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1c.fgh<-function(x) {
   f<-misra1c.f(x)
   g<-misra1c.g(x)
   H<-misra1c.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1<-c(500, 0.0001)
names(start1) <- c("b1", "b2")

start2<-c(600, 0.0002)
names(start2) <- c("b1", "b2")

## Examples

library(numDeriv)

cat("f - start1 =", misra1c.f(start1), "\n")
JA <- misra1c.jac(start1)
JN <- jacobian(misra1c.res, start1)
cat("max(abs(diff)) =", max(abs(JA-JN),"\n"))

cat("f - start2 =", misra1c.f(start2), "\n")
JA <- misra1c.jac(start2)
JN <- jacobian(misra1c.res, start2)
cat("max(abs(diff)) =", max(abs(JA-JN),"\n"))



Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Misra1c)
Try(fm1 <- nls(y ~ b1*(1-(1+2*b2*x)**(-.5)), data = Misra1c, trace = TRUE,
               start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1*(1-(1+2*b2*x)**(-.5)), data = Misra1c, trace = TRUE,
                alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1*(1-(1+2*b2*x)**(-.5)), data = Misra1c, trace = TRUE,
               start = c(b1 = 600, b2 = 0.0002) ))
Try(fm2a <- nls(y ~ b1*(1-(1+2*b2*x)**(-.5)), data = Misra1c, trace = TRUE,
                alg = "port", start = c(b1 = 600, b2 = 0.0002) ))
Try(fm3 <- nls(y ~ 1-(1+2*b2*x)**(-.5), data = Misra1c, trace = TRUE,
               start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ 1-(1+2*b2*x)**(-.5), data = Misra1c, trace = TRUE,
               start = c(b2 = 0.0002), algorithm = "plinear" ))

library(nlsr)

Misra1cnlsr1 <- nlxb(start=start1, formula=Misra1c.formula, data=mypdata, trace=TRUE)
print(Misra1cnlsr1)

Misra1cnlsr2 <- nlxb(start=start2, formula=Misra1c.formula, data=mypdata, trace=TRUE)
print(Misra1cnlsr2)

library(optimr)

Misra1copm1fwd <- opm(start1, misra1c.f, "grfwd", method="ALL")
summary(Misra1copm1fwd, order=value)

Misra1copm2fwd <- opm(start2, misra1c.f, "grfwd", method="ALL")
summary(Misra1copm2fwd, order=value)

Misra1copm1 <- opm(start1, misra1c.f, misra1c.g, method="ALL")
summary(Misra1copm1, order=value)

Misra1copm2 <- opm(start2, misra1c.f, misra1c.g, method="ALL")
summary(Misra1copm2, order=value)

