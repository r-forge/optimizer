## @knitr ##Misra1b.prb
# This is file ##Misra1b.prb
rm(list=ls())
probname <- "##Misra1b"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Misra1b           (Misra1b.dat)

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
Lower Level of Difficulty
Observed Data

Model:         Miscellaneous Class
2 Parameters (b1 and b2)

y = b1 * (1-(1+b2*x/2)**(-2))  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   500         300           3.3799746163E+02  3.1643950207E+00
b2 =     0.0001      0.0002      3.9039091287E-04  4.2547321834E-06

Residual Sum of Squares:                    7.5464681533E-02
Residual Standard Deviation:                7.9301471998E-02
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
Misra1b.formula <- ( y ~ b1 * (1-(1+b2*x/2)**(-2)) )
  
#- setup
  
library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Misra1b")))


misra1b.f <- function(x) {
   res<-misra1b.res(x)
   f<-sum(res*res)
}

misra1b.res <- function(b) {
   xx<-Misra1b$x
   yy<-Misra1b$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]

   res<-b1 * (1-(1+b2*xx/2)**(-2)) - yy
   return(res)
}

# misra1b - Jacobian
misra1b.jac <- function(b) {
stop("not defined")
   xx<-misra1b$x
   yy<-misra1b$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

misra1b.h <- function(x) {
stop("not defined")
   JJ<-misra1b.jac(x)
   H <- t(JJ) %*% JJ
   res<-misra1b.res(x)
stop("not defined")

}

misra1b.g<-function(x) {
#   stop("not defined")
   JJ<-misra1b.jac(x)
   res<-misra1b.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

misra1b.fgh<-function(x) {
   f<-misra1b.f(x)
   g<-misra1b.g(x)
   H<-misra1b.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}


start1<-c(500, 0.0001)
names(start1) <- c("b1", "b2")
start2<-c(300, 0.0002)
names(start2) <- c("b1", "b2")

## Examples
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Misra1b)
Try(fm1 <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
               start = c(b1 = 500, b2 = 0.0001) ))
Try(fm1a <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
                alg = "port", start = c(b1 = 500, b2 = 0.0001) ))
Try(fm2 <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
               start = c(b1 = 300, b2 = 0.0002) ))
Try(fm2a <- nls(y ~ b1 * (1-(1+b2*x/2)**(-2)), data = Misra1b, trace = TRUE,
                alg = "port", start = c(b1 = 300, b2 = 0.0002) ))
Try(fm3 <- nls(y ~ 1-(1+b2*x/2)**(-2), data = Misra1b, trace = TRUE,
               start = c(b2 = 0.0001), algorithm = "plinear" ))
Try(fm4 <- nls(y ~ 1-(1+b2*x/2)**(-2), data = Misra1b, trace = TRUE,
               start = c(b2 = 0.0005), algorithm = "plinear" ))

library(nlsr)

Misra1bnlsr1 <- nlxb(start=start1, formula=Misra1b.formula, data=mypdata, trace=TRUE)
print(Misra1bnlsr1)

Misra1bnlsr2 <- nlxb(start=start2, formula=Misra1b.formula, data=mypdata, trace=TRUE)
print(Misra1bnlsr2)
