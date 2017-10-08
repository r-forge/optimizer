## @knitr ##MGH10.prb
# This is file ##MGH10.prb
rm(list=ls())
probname <- "##MGH10"
probdesc <- "
NIST/ITL StRD
Dataset Name:  MGH10             (MGH10.dat)

File Format:   ASCII
Starting Values   (lines 41 to 43)
Certified Values  (lines 41 to 48)
Data              (lines 61 to 76)

Procedure:     Nonlinear Least Squares Regression

Description:   This problem was found to be difficult for some very
good algorithms.

See More, J. J., Garbow, B. S., and Hillstrom, K. E. 
(1981).  Testing unconstrained optimization software.
ACM Transactions on Mathematical Software. 7(1): 
pp. 17-41.

Reference:     Meyer, R. R. (1970).  
Theoretical and computational aspects of nonlinear 
regression.  In Nonlinear Programming, Rosen, 
Mangasarian and Ritter (Eds).  
New York, NY: Academic Press, pp. 465-486.

Data:          1 Response  (y)
1 Predictor (x)
16 Observations
Higher Level of Difficulty
Generated Data

Model:         Exponential Class
3 Parameters (b1 to b3)

y = b1 * exp[b2/(x+b3)]  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =        2         0.02       5.6096364710E-03  1.5687892471E-04
b2 =   400000      4000          6.1813463463E+03  2.3309021107E+01
b3 =    25000       250          3.4522363462E+02  7.8486103508E-01

Residual Sum of Squares:                    8.7945855171E+01
Residual Standard Deviation:                2.6009740065E+00
Degrees of Freedom:                                13
Number of Observations:                            16











Data:  y               x
3.478000E+04    5.000000E+01
2.861000E+04    5.500000E+01
2.365000E+04    6.000000E+01
1.963000E+04    6.500000E+01
1.637000E+04    7.000000E+01
1.372000E+04    7.500000E+01
1.154000E+04    8.000000E+01
9.744000E+03    8.500000E+01
8.261000E+03    9.000000E+01
7.030000E+03    9.500000E+01
6.005000E+03    1.000000E+02
5.147000E+03    1.050000E+02
4.427000E+03    1.100000E+02
3.820000E+03    1.150000E+02
3.307000E+03    1.200000E+02
2.872000E+03    1.250000E+02

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
MGH10.formula <- ( y ~ b1 * exp(b2/(x+b3)) )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("MGH10")))
# Optimization test function mgh10
# mgh10 from NISTnls
# ??ref...


mgh10.f <- function(x) {
   res<-mgh10.res(x)
   f<-sum(res*res)
}

mgh10.res <- function(b) {
   xx<-MGH10$x # note case!
   yy<-MGH10$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<-b1 * exp(b2/(xx+b3)) - yy
   return(res)
}

# mgh10 - Jacobian
mgh10.jac <- function(b) {
stop("not defined")
   xx<-mgh10$x
   yy<-mgh10$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

mgh10.h <- function(x) {
stop("not defined")
   JJ<-mgh10.jac(x)
   H <- t(JJ) %*% JJ
   res<-mgh10.res(x)
stop("not defined")

}

mgh10.g<-function(x) {
#   stop("not defined")
   JJ<-mgh10.jac(x)
   res<-mgh10.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

mgh10.fgh<-function(x) {
   f<-mgh10.f(x)
   g<-mgh10.g(x)
   H<-mgh10.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}



start1 <- c(b1 = 2, b2 = 400000, b3 = 25000)
start2 <- c(b1 = 0.02, b2 = 4000, b3 = 250)


Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = MGH10)
## check plot on log scale for shape
plot(y ~ x, data = MGH10, log = "y")
## starting values for this run are ridiculous
Try(fm1 <- nls(y ~ b1 * exp(b2/(x+b3)), data = MGH10, trace = TRUE,
               start = c(b1 = 2, b2 = 400000, b3 = 25000)))
Try(fm1a <- nls(y ~ b1 * exp(b2/(x+b3)), data = MGH10,
                trace = TRUE, alg = "port",
                start = c(b1 = 2, b2 = 400000, b3 = 25000)))
Try(fm2 <- nls(y ~ b1 * exp(b2/(x+b3)), data = MGH10, trace = TRUE,
               start = c(b1 = 0.02, b2 = 4000, b3 = 250)))
Try(fm2a <- nls(y ~ b1 * exp(b2/(x+b3)), data = MGH10,
                trace = TRUE, alg = "port",
                start = c(b1 = 0.02, b2 = 4000, b3 = 250)))
Try(fm3 <- nls(y ~ exp(b2/(x+b3)), data = MGH10, trace = TRUE,
               start = c(b2 = 400000, b3 = 25000),
               algorithm = "plinear"))
Try(fm4 <- nls(y ~ exp(b2/(x+b3)), data = MGH10, trace = TRUE,
               start = c(b2 = 4000, b3 = 250),
               algorithm = "plinear"))

library(nlsr)

MGH10nlxb1 <- nlxb(start=start1, formula=MGH10.formula, data=mypdata, trace=TRUE)
print(MGH10nlxb1)

MGH10nlxb2 <- nlxb(start=start2, formula=MGH10.formula, data=mypdata, trace=TRUE)
print(MGH10nlxb2)

library(optimr)
mset <- c("Rvmmin", "bobyqa", "L-BFGS-B")
MGH10opm1fwd <- opm(start1, mgh10.f, "grfwd", method=mset)
summary(MGH10opm1fwd, order=value)

MGH10opm2fwd <- opm(start2, mgh10.f, "grfwd", method=mset)
summary(MGH10opm2fwd, order=value)

MGH10opm1central <- opm(start1, mgh10.f, "grcentral", method=mset)
summary(MGH10opm1central, order=value)

MGH10opm2central <- opm(start2, mgh10.f, "grcentral", method=mset)
summary(MGH10opm2central, order=value)

