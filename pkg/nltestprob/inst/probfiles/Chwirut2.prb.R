## @knitr ##Chwirut2.prb
# This is file ##Chwirut2.prb
probname <- "##Chwirut2"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Chwirut2          (Chwirut2.dat)

File Format:   ASCII
Starting Values   (lines 41 to  43)
Certified Values  (lines 41 to  48)
Data              (lines 61 to 114)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study involving
ultrasonic calibration.  The response variable is
ultrasonic response, and the predictor variable is
metal distance.



Reference:     Chwirut, D., NIST (197?).  
Ultrasonic Reference Block Study. 





Data:          1 Response  (y = ultrasonic response)
1 Predictor (x = metal distance)
54 Observations
Lower Level of Difficulty
Observed Data

Model:         Exponential Class
3 Parameters (b1 to b3)

y = exp(-b1*x)/(b2+b3*x)  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   0.1         0.15          1.6657666537E-01  3.8303286810E-02
b2 =   0.01        0.008         5.1653291286E-03  6.6621605126E-04
b3 =   0.02        0.010         1.2150007096E-02  1.5304234767E-03

Residual Sum of Squares:                    5.1304802941E+02
Residual Standard Deviation:                3.1717133040E+00
Degrees of Freedom:                                51
Number of Observations:                            54











Data:  y             x
92.9000E0     0.500E0
57.1000E0     1.000E0
31.0500E0     1.750E0
11.5875E0     3.750E0
8.0250E0     5.750E0
63.6000E0     0.875E0
21.4000E0     2.250E0
14.2500E0     3.250E0
8.4750E0     5.250E0
63.8000E0     0.750E0
26.8000E0     1.750E0
16.4625E0     2.750E0
7.1250E0     4.750E0
67.3000E0     0.625E0
41.0000E0     1.250E0
21.1500E0     2.250E0
8.1750E0     4.250E0
81.5000E0      .500E0
13.1200E0     3.000E0
59.9000E0      .750E0
14.6200E0     3.000E0
32.9000E0     1.500E0
5.4400E0     6.000E0
12.5600E0     3.000E0
5.4400E0     6.000E0
32.0000E0     1.500E0
13.9500E0     3.000E0
75.8000E0      .500E0
20.0000E0     2.000E0
10.4200E0     4.000E0
59.5000E0      .750E0
21.6700E0     2.000E0
8.5500E0     5.000E0
62.0000E0      .750E0
20.2000E0     2.250E0
7.7600E0     3.750E0
3.7500E0     5.750E0
11.8100E0     3.000E0
54.7000E0      .750E0
23.7000E0     2.500E0
11.5500E0     4.000E0
61.3000E0      .750E0
17.7000E0     2.500E0
8.7400E0     4.000E0
59.2000E0      .750E0
16.3000E0     2.500E0
8.6200E0     4.000E0
81.0000E0      .500E0
4.8700E0     6.000E0
14.6200E0     3.000E0
81.7000E0      .500E0
17.1700E0     2.750E0
81.3000E0      .500E0
28.9000E0     1.750E0

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
##Chwirut2.formula <- ( y ~ b1*x**b2 )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Chwirut2")))

# Optimization test function chwirut2
# chwirut2 from NISTnls
# ??ref...


chwirut2.f <- function(x) {
   res<-chwirut2.res(x)
   f<-sum(res*res)
}

chwirut2.res <- function(b) {
   xx<-Chwirut2$x # note case !!
   yy<-Chwirut2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-exp(-b1*xx)/(b2+b3*xx) - yy
   return(res)
}

# chwirut2 - Jacobian
chwirut2.jac <- function(b) {
xx<-Chwirut2$x # note caps!
J <- matrix(0, length(xx), 3) # define the size of the Jacobian
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
#   res<-b1*(b2+xx)**(-1/b3) - yy
   exx<-exp(-b1*xx)
   den<-(b2+b3*xx)
   expr1 <- b2 + xx
    expr3 <- -1/b3
    expr4 <- expr1^expr3
    J[ , 1] <- -xx*exx/den
    J[ , 2] <- -exx/(den*den)
    J[ , 3] <- -exx*xx/(den*den)
return(J)
}

chwirut2.h <- function(x) {
stop("not defined")
   JJ<-chwirut2.jac(x)
   H <- t(JJ) %*% JJ
   res<-chwirut2.res(x)
stop("not defined")

}

chwirut2.g<-function(x) {
#   stop("not defined")
   JJ<-chwirut2.jac(x)
   res<-chwirut2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

chwirut2.fgh<-function(x) {
   f<-chwirut2.f(x)
   g<-chwirut2.g(x)
   H<-chwirut2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

chwirut2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Chwirut2) # and load up the data into x and y
   start<-c( 0.1, 0.01, 0.02)
}

chwirut2.test<-function() {
}   

## Examples

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Chwirut2)
Try(fm1 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02)))
Try(fm1a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
                start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02), alg = "port"))
Try(fm2 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01)))
Try(fm2a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
                start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01), alg = "port"))
Try(fm3 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.1, p3 = 2.), alg = "plinear"))
Try(fm4 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.15, p3 = 0.01/0.008), alg = "plinear"))
