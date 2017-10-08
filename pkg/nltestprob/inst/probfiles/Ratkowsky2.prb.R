## @knitr ##Ratkowsky2
# This is file ##Ratkowsky2.prb.R
rm(list=ls())
probname <- "##Ratkowsky2"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Rat42             (Rat42.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 43)
               Certified Values  (lines 41 to 48)
               Data              (lines 61 to 69)

Procedure:     Nonlinear Least Squares Regression

Description:   This model and data are an example of fitting
               sigmoidal growth curves taken from Ratkowsky (1983).
               The response variable is pasture yield, and the
               predictor variable is growing time.


Reference:     Ratkowsky, D.A. (1983).  
               Nonlinear Regression Modeling.
               New York, NY:  Marcel Dekker, pp. 61 and 88.





Data:          1 Response  (y = pasture yield)
               1 Predictor (x = growing time)
               9 Observations
               Higher Level of Difficulty
               Observed Data

Model:         Exponential Class
               3 Parameters (b1 to b3)

               y = b1 / (1+exp[b2-b3*x])  +  e



          Starting Values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   100         75            7.2462237576E+01  1.7340283401E+00
  b2 =     1          2.5          2.6180768402E+00  8.8295217536E-02
  b3 =     0.1        0.07         6.7359200066E-02  3.4465663377E-03

Residual Sum of Squares:                    8.0565229338E+00
Residual Standard Deviation:                1.1587725499E+00
Degrees of Freedom:                                6
Number of Observations:                            9 











Data:   y              x
       8.930E0        9.000E0
      10.800E0       14.000E0
      18.590E0       21.000E0
      22.330E0       28.000E0
      39.350E0       42.000E0
      56.110E0       57.000E0
      61.730E0       63.000E0
      64.620E0       70.000E0
      67.080E0       79.000E0

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
ratkowsky2.formula <- ( y ~ b1 / (1+exp(b2-b3*x)) )


ratkowsky2.f <- function(x) {
   res<-ratkowsky2.res(x)
   f<-sum(res*res)
}

ratkowsky2.res <- function(b) {
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-b1 / (1+exp(b2-b3*xx)) - yy
   return(res)
}

# ratkowsky2 - Jacobian
ratkowsky2.jac <- function(b) {
stop("not defined")
   xx<-Ratkowsky2$x
   yy<-Ratkowsky2$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

ratkowsky2.h <- function(x) {
stop("not defined")
   JJ<-ratkowsky2.jac(x)
   H <- t(JJ) %*% JJ
   res<-ratkowsky2.res(x)
stop("not defined")

}

ratkowsky2.g<-function(x) {
#   stop("not defined")
   JJ<-ratkowsky2.jac(x)
   res<-ratkowsky2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

ratkowsky2.fgh<-function(x) {
   f<-ratkowsky2.f(x)
   g<-ratkowsky2.g(x)
   H<-ratkowsky2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}


   start1<-c(100,1,.1)
   start2<-c(75,2.5,0.07)
   start0<-rep(1,3)
   names(start0) <- c("b1","b2","b3")
   names(start1) <- c("b1","b2","b3")
   names(start2) <- c("b1","b2","b3")
   
library("NISTnls", character.only=TRUE) # get parent collection

mypdata <- eval(parse(text=data(list="Ratkowsky2")))
cat("Rat42=Ratkowsky2 data:\n")
print(mypdata)

cat("nls tries\n")
rat2nls0 <- try(nls(formula=ratkowsky2.formula, start=start0, trace=TRUE, data=mypdata))
print(rat2nls0)

rat2nls1 <- try(nls(formula=ratkowsky2.formula, start=start1, trace=TRUE, data=mypdata))
print(rat2nls1)

rat2nls2 <- try(nls(formula=ratkowsky2.formula, start=start2, trace=TRUE, data=mypdata))
print(rat2nls2)

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Ratkowsky2)

Try(fm1 <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
               start = c(b1 = 100, b2 = 1, b3 = 0.1)))
Try(fm1a <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2,
                trace = TRUE, alg = "port", 
                start = c(b1 = 100, b2 = 1, b3 = 0.1)))
Try(fm2 <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
               start = c(b1 = 75, b2 = 2.5, b3 = 0.07)))
Try(fm2a <- nls(y ~ b1 / (1+exp(b2-b3*x)), data = Ratkowsky2,
                trace = TRUE, alg = "port", 
                start = c(b1 = 75, b2 = 2.5, b3 = 0.07)))
Try(fm3 <- nls(y ~ 1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
               start = c(b2 = 1, b3 = 0.1), alg = "plinear"))
Try(fm4 <- nls(y ~ 1 / (1+exp(b2-b3*x)), data = Ratkowsky2, trace = TRUE,
               start = c(b2 = 2.5, b3 = 0.07), alg = "plinear"))

## Using a self-starting model
Try(fm5 <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = Ratkowsky2))
summary(fm5)

cat("nlxb tries\n")
library(nlsr)
rat2nlxb0 <- try(nlxb(formula=ratkowsky2.formula, start=start0, trace=TRUE, data=mypdata))
print(rat2nlxb0)

rat2nlxb1 <- try(nlxb(formula=ratkowsky2.formula, start=start1, trace=TRUE, data=mypdata))
print(rat2nlxb1)

rat2nlxb2 <- try(nlxb(formula=ratkowsky2.formula, start=start2, trace=TRUE, data=mypdata))
print(rat2nlxb2)

cat("optimr tries\n")
library(optimr)
rat2opm0 <- opm(start0, ratkowsky2.f, "grcentral", method="ALL")
summary(rat2opm0, order=value)

rat2opm1 <- opm(start1, ratkowsky2.f, "grcentral", method="ALL")
summary(rat2opm1, order=value)

rat2opm2 <- opm(start2, ratkowsky2.f, "grcentral", method="ALL")
summary(rat2opm2, order=value)

