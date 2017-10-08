## @knitr ##Nelson.prb
# This is file ##Nelson.prb
rm(list=ls())
probname <- "##Nelson"
probdesc <- "

NIST/ITL StRD
Dataset Name:  Nelson            (Nelson.dat)

File Format:   ASCII
Starting Values   (lines 41 to 43)
Certified Values  (lines 41 to 48)
Data              (lines 61 to 188)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a study involving
the analysis of performance degradation data from
accelerated tests, published in IEEE Transactions
on Reliability.  The response variable is dialectric
breakdown strength in kilo-volts, and the predictor
variables are time in weeks and temperature in degrees
Celcius.


Reference:     Nelson, W. (1981).  
Analysis of Performance-Degradation Data.  
IEEE Transactions on Reliability.
Vol. 2, R-30, No. 2, pp. 149-155.

Data:          1 Response   ( y = dialectric breakdown strength) 
2 Predictors (x1 = time; x2 = temperature)
128 Observations
Average Level of Difficulty
Observed Data

Model:         Exponential Class
3 Parameters (b1 to b3)

log[y] = b1 - b2*x1 * exp[-b3*x2]  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =    2           2.5          2.5906836021E+00  1.9149996413E-02
b2 =    0.0001      0.000000005  5.6177717026E-09  6.1124096540E-09
b3 =   -0.01       -0.05        -5.7701013174E-02  3.9572366543E-03

Residual Sum of Squares:                    3.7976833176E+00
Residual Standard Deviation:                1.7430280130E-01
Degrees of Freedom:                               125
Number of Observations:                           128











Data:   y              x1            x2
15.00E0         1E0         180E0
17.00E0         1E0         180E0
15.50E0         1E0         180E0
16.50E0         1E0         180E0
15.50E0         1E0         225E0
15.00E0         1E0         225E0
16.00E0         1E0         225E0
14.50E0         1E0         225E0
15.00E0         1E0         250E0
14.50E0         1E0         250E0
12.50E0         1E0         250E0
11.00E0         1E0         250E0
14.00E0         1E0         275E0
13.00E0         1E0         275E0
14.00E0         1E0         275E0
11.50E0         1E0         275E0
14.00E0         2E0         180E0
16.00E0         2E0         180E0
13.00E0         2E0         180E0
13.50E0         2E0         180E0
13.00E0         2E0         225E0
13.50E0         2E0         225E0
12.50E0         2E0         225E0
12.50E0         2E0         225E0
12.50E0         2E0         250E0
12.00E0         2E0         250E0
11.50E0         2E0         250E0
12.00E0         2E0         250E0
13.00E0         2E0         275E0
11.50E0         2E0         275E0
13.00E0         2E0         275E0
12.50E0         2E0         275E0
13.50E0         4E0         180E0
17.50E0         4E0         180E0
17.50E0         4E0         180E0
13.50E0         4E0         180E0
12.50E0         4E0         225E0
12.50E0         4E0         225E0
15.00E0         4E0         225E0
13.00E0         4E0         225E0
12.00E0         4E0         250E0
13.00E0         4E0         250E0
12.00E0         4E0         250E0
13.50E0         4E0         250E0
10.00E0         4E0         275E0
11.50E0         4E0         275E0
11.00E0         4E0         275E0
9.50E0         4E0         275E0
15.00E0         8E0         180E0
15.00E0         8E0         180E0
15.50E0         8E0         180E0
16.00E0         8E0         180E0
13.00E0         8E0         225E0
10.50E0         8E0         225E0
13.50E0         8E0         225E0
14.00E0         8E0         225E0
12.50E0         8E0         250E0
12.00E0         8E0         250E0
11.50E0         8E0         250E0
11.50E0         8E0         250E0
6.50E0         8E0         275E0
5.50E0         8E0         275E0
6.00E0         8E0         275E0
6.00E0         8E0         275E0
18.50E0        16E0         180E0
17.00E0        16E0         180E0
15.30E0        16E0         180E0
16.00E0        16E0         180E0
13.00E0        16E0         225E0
14.00E0        16E0         225E0
12.50E0        16E0         225E0
11.00E0        16E0         225E0
12.00E0        16E0         250E0
12.00E0        16E0         250E0
11.50E0        16E0         250E0
12.00E0        16E0         250E0
6.00E0        16E0         275E0
6.00E0        16E0         275E0
5.00E0        16E0         275E0
5.50E0        16E0         275E0
12.50E0        32E0         180E0
13.00E0        32E0         180E0
16.00E0        32E0         180E0
12.00E0        32E0         180E0
11.00E0        32E0         225E0
9.50E0        32E0         225E0
11.00E0        32E0         225E0
11.00E0        32E0         225E0
11.00E0        32E0         250E0
10.00E0        32E0         250E0
10.50E0        32E0         250E0
10.50E0        32E0         250E0
2.70E0        32E0         275E0
2.70E0        32E0         275E0
2.50E0        32E0         275E0
2.40E0        32E0         275E0
13.00E0        48E0         180E0
13.50E0        48E0         180E0
16.50E0        48E0         180E0
13.60E0        48E0         180E0
11.50E0        48E0         225E0
10.50E0        48E0         225E0
13.50E0        48E0         225E0
12.00E0        48E0         225E0
7.00E0        48E0         250E0
6.90E0        48E0         250E0
8.80E0        48E0         250E0
7.90E0        48E0         250E0
1.20E0        48E0         275E0
1.50E0        48E0         275E0
1.00E0        48E0         275E0
1.50E0        48E0         275E0
13.00E0        64E0         180E0
12.50E0        64E0         180E0
16.50E0        64E0         180E0
16.00E0        64E0         180E0
11.00E0        64E0         225E0
11.50E0        64E0         225E0
10.50E0        64E0         225E0
10.00E0        64E0         225E0
7.27E0        64E0         250E0
7.50E0        64E0         250E0
6.70E0        64E0         250E0
7.60E0        64E0         250E0
1.50E0        64E0         275E0
1.00E0        64E0         275E0
1.20E0        64E0         275E0
1.20E0        64E0         275E0
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

nelson.formula <- ( log(y) ~ b1 - b2*x1 * exp(-b3*x2) )


nelson.f <- function(x) {
   res<-nelson.res(x)
   f<-sum(res*res)
}

nelson.res <- function(b) {
   xx1<-mypdata$x1
   xx2<-mypdata$x2
   yy<-mypdata$y
   logyy<-log(yy)
   res <- rep(NA, length(yy))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]

   res<- b1 - b2*xx1 * exp(-b3*xx2) - logyy
   return(res)
}

# nelson - Jacobian
nelson.jac <- function(b) {
  stop("not defined yet")
   xx1<-mypdata$x1
   xx2<-mypdata$x2
   yy<-mypdata$y
   logyy<-log(yy)
   n<-length(b)
   m<-length(yy)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
   return(J)
}

nelson.h <- function(x) {
stop("not defined")
   JJ<-nelson.jac(x)
   H <- t(JJ) %*% JJ
   res<-nelson.res(x)
stop("not defined")

}

nelson.g<-function(x) {
#   stop("not defined")
   JJ<-nelson.jac(x)
   res<-nelson.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

nelson.fgh<-function(x) {
   f<-nelson.f(x)
   g<-nelson.g(x)
   H<-nelson.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

library(NISTnls) # get parent collection
mypdata <- eval(parse(text=data(list="Nelson")))
cat("Nelson problem:\n")
print(mypdata)


start1 <- c(2,0.0001,-0.01)
start2 <- c(2.5, 0.000000005, -0.05)
start0 <- rep(0,3)
names(start0) <- c("b1", "b2", "b3")
names(start1) <- names(start0)
names(start2) <- names(start0)

nelson0nls <- try(nls(nelson.formula, start=start0, trace=TRUE, data=mypdata))
print(nelson0nls)

nelson1nls <- try(nls(nelson.formula, start=start1, trace=TRUE, data=mypdata))
print(nelson1nls)

nelson2nls <- try(nls(nelson.formula, start=start2, trace=TRUE, data=mypdata))
print(nelson2nls)


Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x1, data = Nelson, log = "y")
plot(y ~ x2, data = Nelson, log = "y")
coplot(y ~ x1 | x2, data = Nelson)
coplot(y ~ x2 | x1, data = Nelson)

Try(fm1 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
               start = c(b1 = 2, b2 = 0.0001, b3 = -0.01), trace = TRUE))
Try(fm1a <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
                trace = TRUE, alg = "port",
                start = c(b1 = 2, b2 = 0.0001, b3 = -0.01)))

Try(fm2 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
               start = c(b1 = 2.5, b2 = 0.000000005, b3 = -0.05), trace = TRUE))
Try(fm2 <- nls(log(y) ~ b1 - b2*x1 * exp(-b3*x2), data = Nelson,
               trace = TRUE, alg = "port", 
               start = c(b1 = 2.5, b2 = 0.000000005, b3 = -0.05)))

Try(fm3 <- nls(log(y) ~ cbind(1, -x1 * exp(-b3*x2)), data = Nelson,
               start = c(b3 = -0.01), trace = TRUE, algorithm = "plinear"))

Try(fm4 <- nls(log(y) ~ cbind(1, -x1 * exp(-b3*x2)), data = Nelson,
               start = c(b3 = -0.05), trace = TRUE, algorithm = "plinear"))


nelson0nlxb <- nlxb(nelson.formula, start=start0, trace=TRUE, data=mypdata)
print(nelson0nlxb)

nelson1nlxb <- nlxb(nelson.formula, start=start1, trace=TRUE, data=mypdata)
print(nelson1nlxb)

nelson2nlxb <- nlxb(nelson.formula, start=start2, trace=TRUE, data=mypdata)
print(nelson2nlxb)

cat("Test Nelson R functions\n")
cat("start0 =")
print(start0)
cat("function value at start0=", nelson.f(start0))
cat("start1 =")
print(start1)
cat("function value at start1=", nelson.f(start1))
cat("start2 =")
print(start2)
cat("function value at start2=", nelson.f(start2))

tobedefined <- '

library(numDeriv)

ga0 <- nelson.g(start0)
gn0 <- grad(nelson.f, start0)
cat("Gradient at start0 =")
print(ga0)
cat("Numerical approximation at start0 =")
print(gn0)
cat("max(abs(diff))= ", max(abs(ga0-gn0)), "\n")


ga1 <- nelson.g(start1)
gn1 <- grad(nelson.f, start1)
cat("Gradient at start1 =")
print(ga1)
cat("Numerical approximation at start1 =")
print(gn1)
cat("max(abs(diff))= ", max(abs(ga1-gn1)), "\n")

ga2 <- nelson.g(start2)
gn2 <- grad(nelson.f, start2)
cat("Gradient at start2 =")
print(ga2)
cat("Numerical approximation at start2 =")
print(gn2)
cat("max(abs(diff))= ", max(abs(ga2-gn2)), "\n")
'
## REMOVE QUOTES WHEN FIXED

library(optimr)

nelson0opmnumg <- opm(start0, nelson.f, method="ALL")
summary(nelson0opmnumg, order=value)
