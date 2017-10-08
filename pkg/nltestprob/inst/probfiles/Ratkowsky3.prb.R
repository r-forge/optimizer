## @knitr ##Ratkowsky3.prb
# This is file ##Ratkowsky3.prb.R
rm(list=ls())

probname <- "##Ratkowsky3"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Ratkowsky3        (Ratkowsky3.dat)

File Format:   ASCII
Starting Values   (lines 41 to 44)
Certified Values  (lines 41 to 49)
Data              (lines 61 to 75)

Procedure:     Nonlinear Least Squares Regression

Description:   This model and data are an example of fitting  
sigmoidal growth curves taken from Ratkowsky (1983).  
The response variable is the dry weight of onion bulbs 
and tops, and the predictor variable is growing time. 


Reference:     Ratkowsky, D.A. (1983).  
Nonlinear Regression Modeling.
New York, NY:  Marcel Dekker, pp. 62 and 88.





Data:          1 Response  (y = onion bulb dry weight)
1 Predictor (x = growing time)
15 Observations
Higher Level of Difficulty
Observed Data

Model:         Exponential Class
4 Parameters (b1 to b4)

y = b1 / ((1+exp[b2-b3*x])**(1/b4))  +  e



Starting Values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   100         700           6.9964151270E+02  1.6302297817E+01
b2 =    10           5           5.2771253025E+00  2.0828735829E+00
b3 =     1           0.75        7.5962938329E-01  1.9566123451E-01
b4 =     1           1.3         1.2792483859E+00  6.8761936385E-01

Residual Sum of Squares:                    8.7864049080E+03
Residual Standard Deviation:                2.8262414662E+01
Degrees of Freedom:                                9
Number of Observations:                           15 










Data:   y          x
16.08E0     1.0E0
33.83E0     2.0E0
65.80E0     3.0E0
97.20E0     4.0E0
191.55E0     5.0E0
326.20E0     6.0E0
386.87E0     7.0E0
520.53E0     8.0E0
590.03E0     9.0E0
651.92E0    10.0E0
724.93E0    11.0E0
699.56E0    12.0E0
689.96E0    13.0E0
637.56E0    14.0E0
717.41E0    15.0E0
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
ratkowsky3.formula <- ( y ~  b1 / ((1+exp(b2-b3*x))**(1/b4)) ) 
# Optimization test function ratkowsky3
# ratkowsky3 from NISTnls
# ??ref...

library("NISTnls", character.only=TRUE) # get parent collection

mypdata <- eval(parse(text=data(list="Ratkowsky3")))



ratkowsky3.f <- function(x) {
   res<-ratkowsky3.res(x)
   f<-sum(res*res)
}

ratkowsky3.res <- function(b) {
   xx<-Ratkowsky3$x
   yy<-Ratkowsky3$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   res<- b1 / ((1+exp(b2-b3*xx))**(1/b4)) - yy
   return(res)
}

# ratkowsky3 - Jacobian
ratkowsky3.jac <- function(b) {
stop("not defined")
   xx<-Ratkowsky3$x
   yy<-Ratkowsky3$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

ratkowsky3.h <- function(x) {
stop("not defined")
   JJ<-ratkowsky3.jac(x)
   H <- t(JJ) %*% JJ
   res<-ratkowsky3.res(x)
stop("not defined")

}

ratkowsky3.g<-function(x) {
#   stop("not defined")
   JJ<-ratkowsky3.jac(x)
   res<-ratkowsky3.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

ratkowsky3.fgh<-function(x) {
   f<-ratkowsky3.f(x)
   g<-ratkowsky3.g(x)
   H<-ratkowsky3.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

start1<-c(100,10,1,1)
names(start1) <- c("b1","b2","b3","b4")
start2<-c(700,5,0.75,1.3)
names(start2) <- c("b1","b2","b3","b4")
start0<-rep(1,4)
names(start0) <- c("b1","b2","b3","b4")

NISTRat3<-list()
NISTRat3$value<-8.7864049080E+03
NISTRat3$par<-c(6.9964151270E+02,5.2771253025E+00,7.5962938329E-01,1.2792483859E+00)
NISTRat3$ses<-c(1.6302297817E+01, 2.0828735829E+00, 1.9566123451E-01,6.8761936385E-01)



Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
## causes NA/NaN/Inf error
Try(fm1 <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
               start = c(b1 = 100, b2 = 10, b3 = 1, b4 = 1),
               trace = TRUE))
Try(fm1a <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
                start = c(b1 = 100, b2 = 10, b3 = 1, b4 = 1),
                alg = "port", trace = TRUE))

Try(fm2 <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
               start = c(b1 = 700, b2 = 5, b3 = 0.75, b4 = 1.3),
               trace = TRUE))
Try(fm2a <- nls(y ~ b1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
                start = c(b1 = 700, b2 = 5, b3 = 0.75, b4 = 1.3),
                alg = "port", trace = TRUE))

Try(fm3 <- nls(y ~ 1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
               start = c(b2 = 10, b3 = 1, b4 = 1), algorithm = "plinear",
               trace = TRUE))
Try(fm4 <- nls(y ~ 1 / ((1+exp(b2-b3*x))**(1/b4)), data = Ratkowsky3,
               start = c(b2 = 5, b3 = 0.75, b4 = 1.3), algorithm = "plinear",
               trace = TRUE))

# nlxb

library(nlsr)

Rat3nlxb0 <- nlxb(start=start0, formula=ratkowsky3.formula, trace=TRUE, data=mypdata)
print(Rat3nlxb0)

Rat3nlxb1 <- nlxb(start=start1, formula=ratkowsky3.formula, trace=TRUE, data=mypdata)
print(Rat3nlxb1)

Rat3nlxb2 <- nlxb(start=start2, formula=ratkowsky3.formula, trace=TRUE, data=mypdata)
print(Rat3nlxb2)

