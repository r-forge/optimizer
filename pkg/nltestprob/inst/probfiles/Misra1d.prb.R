## @knitr ##Misra1d.prb
# This is file ##Misra1d.prb
rm(list=ls())
probname <- "##Misra1d"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Misra1d           (Misra1d.dat)

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







Data:          1 Response  (y = volume)
1 Predictor (x = pressure)
14 Observations
Average Level of Difficulty
Observed Data

Model:         Miscellaneous Class
2 Parameters (b1 and b2)

y = b1*b2*x*((1+b2*x)**(-1))  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   500         450           4.3736970754E+02  3.6489174345E+00
b2 =     0.0001      0.0003      3.0227324449E-04  2.9334354479E-06

Residual Sum of Squares:                    5.6419295283E-02
Residual Standard Deviation:                6.8568272111E-02
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
Misra1d.formula <- ( y ~ b1*b2*x*((1+b2*x)**(-1)) )

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Misra1c")))
# Optimization test function misra1c

# Optimization test function misra1d
# misra1d from NISTnls
# ??ref...


misra1d.f <- function(x) {
   res<-misra1d.res(x)
   f<-sum(res*res)
}

misra1d.res <- function(b) {
   xx<-Misra1d$x
   yy<-Misra1d$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1*b2*xx*((1+b2*xx)**(-1)) - yy
   return(res)
}

# misra1d - Jacobian

misra1d.jac <- function(b) {
   xx<-Misra1d$x
   yy<-Misra1d$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   J<-matrix(0,m,n) # define the size of the Jacobian
   J[,1] <- b2*xx*((1+b2*xx)**(-1)) 
   J[,2] <- b1* (xx*((1+b2*xx)**(-1)) - b2*xx*xx*((1+b2*xx)**(-2)) )
   J
}

misra1d.h <- function(x) {
stop("not defined")
   JJ<-Misra1d.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1d.res(x)
stop("not defined")

}

misra1d.g<-function(x) {
#   stop("not defined")
   JJ<-misra1d.jac(x)
   res<-misra1d.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1d.fgh<-function(x) {
   f<-misra1d.f(x)
   g<-misra1d.g(x)
   H<-misra1d.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1<-c(500, 0.0001)
names(start1) <- c("b1","b2")
start2<-c(450, 0.0003)
names(start2) <- c("b1","b2")

# Examples

library(numDeriv)

cat("f - start1 =", misra1d.f(start1), "\n")
JA <- misra1d.jac(start1)
JN <- jacobian(misra1d.res, start1)
cat("max(abs(diff)) =", max(abs(JA-JN),"\n"))

cat("f - start2 =", misra1d.f(start2), "\n")
JA <- misra1d.jac(start2)
JN <- jacobian(misra1d.res, start2)
cat("max(abs(diff)) =", max(abs(JA-JN),"\n"))




Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Misra1d)
Try(fm1 <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
               start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
                alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
               start = c(b1 = 450, b2 = 0.0003) ))
Try(fm2a <- nls(y ~ b1*b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
                alg = "port", start = c(b1 = 450, b2 = 0.0003) ))
Try(fm3 <- nls(y ~ b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
               start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ b2*x*((1+b2*x)**(-1)), data = Misra1d, trace = TRUE,
               start = c(b2 = 0.0005), algorithm = "plinear" ))

require(nlsr)

Misra1dnlxb1 <- nlxb(start=start1, formula=Misra1d.formula, data=mypdata, trace=TRUE)
print(Misra1dnlxb1)

Misra1dnlxb2 <- nlxb(start=start2, formula=Misra1d.formula, data=mypdata, trace=TRUE)
print(Misra1dnlxb2)

require(optimr)
   
Misra1dopm1 <- opm(start1, misra1d.f, misra1d.g, method="ALL")
summary(Misra1dopm1, order=value)

Misra1dopm2 <- opm(start2, misra1d.f, misra1d.g, method="ALL")
summary(Misra1dopm2, order=value)
