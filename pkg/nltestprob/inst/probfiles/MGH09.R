## @knitr ##MGH09.prb
# This is file ##MGH09.prb
probname <- "##MGH09"
probdesc <- "NIST/ITL StRD
Dataset Name:  MGH09             (MGH09.dat)

File Format:   ASCII
Starting Values   (lines 41 to 44)
Certified Values  (lines 41 to 49)
Data              (lines 61 to 71)

Procedure:     Nonlinear Least Squares Regression

Description:   This problem was found to be difficult for some very 
good algorithms.  There is a local minimum at (+inf,
-14.07..., -inf, -inf) with final sum of squares 
0.00102734....

See More, J. J., Garbow, B. S., and Hillstrom, K. E. 
(1981).  Testing unconstrained optimization software.
ACM Transactions on Mathematical Software. 7(1): 
pp. 17-41.

Reference:     Kowalik, J.S., and M. R. Osborne, (1978).  
Methods for Unconstrained Optimization Problems.  
New York, NY:  Elsevier North-Holland.

Data:          1 Response  (y)
1 Predictor (x)
11 Observations
Higher Level of Difficulty
Generated Data

Model:         Rational Class (linear/quadratic)
4 Parameters (b1 to b4)

y = b1*(x**2+x*b2) / (x**2+x*b3+b4)  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   25          0.25          1.9280693458E-01  1.1435312227E-02
b2 =   39          0.39          1.9128232873E-01  1.9633220911E-01
b3 =   41.5        0.415         1.2305650693E-01  8.0842031232E-02
b4 =   39          0.39          1.3606233068E-01  9.0025542308E-02

Residual Sum of Squares:                    3.0750560385E-04
Residual Standard Deviation:                6.6279236551E-03
Degrees of Freedom:                                7
Number of Observations:                           11










Data:  y               x
1.957000E-01    4.000000E+00
1.947000E-01    2.000000E+00
1.735000E-01    1.000000E+00
1.600000E-01    5.000000E-01
8.440000E-02    2.500000E-01
6.270000E-02    1.670000E-01
4.560000E-02    1.250000E-01
3.420000E-02    1.000000E-01
3.230000E-02    8.330000E-02
2.350000E-02    7.140000E-02
2.460000E-02    6.250000E-02



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
MGH09.formula <- ( y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4)  )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("MGH09")))


# Optimization test function mgh09
# mgh09 from NISTnls
# ??ref...


mgh09.f <- function(x) {
   res<-mgh09.res(x)
   f<-sum(res*res)
}

mgh09.res <- function(b) {
   xx<-MGH09$x # note case!
   yy<-MGH09$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<-b1*(xx**2+xx*b2) / (xx**2+xx*b3+b4) - yy
   return(res)
}

# mgh09 - Jacobian
mgh09.jac <- function(b) {
stop("not defined")
   xx<-mgh09$x
   yy<-mgh09$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

mgh09.h <- function(x) {
stop("not defined")
   JJ<-mgh09.jac(x)
   H <- t(JJ) %*% JJ
   res<-mgh09.res(x)
stop("not defined")

}

mgh09.g<-function(x) {
#   stop("not defined")
   JJ<-mgh09.jac(x)
   res<-mgh09.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

mgh09.fgh<-function(x) {
   f<-mgh09.f(x)
   g<-mgh09.g(x)
   H<-mgh09.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1 <- c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39) 
start2 = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39)


Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = MGH09)
## starting values for this attempt are ridiculous
Try(fm1 <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
               data = MGH09, trace = TRUE,
               start = c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39)))
Try(fm1a <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
                data = MGH09, trace = TRUE, alg = "port",
                start = c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39)))

Try(fm2 <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
               data = MGH09, trace = TRUE,
               start = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39)))
Try(fm2a <- nls(y ~ b1*(x**2+x*b2) / (x**2+x*b3+b4),
                data = MGH09, trace = TRUE, alg = "port",
                start = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39)))
Try(fm3 <- nls(y ~ cbind(x, x**2) / (x**2+x*b3+b4),
               data = MGH09, trace = TRUE, algorithm = "plinear",
               start = c(b3 = 41.5, b4 = 39)))
Try(fm4 <- nls(y ~ cbind(x, x**2) / (x**2+x*b3+b4),
               data = MGH09, trace = TRUE, algorithm = "plinear",
               start = c(b3 = 0.415, b4 = 0.39)))

# nlsr
library(nlsr)

MGH09nlxb1 <- nlxb(start=start1, formula=MGH09.formula, data=mypdata, trace=TRUE)
print(MGH09nlxb1)
MGH09nlxb2 <- nlxb(start=start2, formula=MGH09.formula, data=mypdata, trace=TRUE)
print(MGH09nlxb2)

library(optimr)

mset <- c("Rvmmin", "bobyqa", "nmkb")

MGH09opm1fwd <- opm(start1, mgh09.f, "grfwd", method=mset)
summary(MGH09opm1fwd, order=value)

MGH09opm2fwd <- opm(start2, mgh09.f, "grfwd", method=mset)
summary(MGH09opm2fwd, order=value)
