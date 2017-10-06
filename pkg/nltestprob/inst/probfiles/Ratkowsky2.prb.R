## @knitr ##Ratkowsky2
# This is file ##Ratkowsky2.prb.R
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
ratkowsky2.formula <- ( y ~ b1*x**b2 )

# Optimization test function ratkowsky2
# ratkowsky2 from NISTnls
# ??ref...


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

ratkowsky2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Ratkowsky2) # and load up the data into x and y
}

ratkowsky2.test<-function() {
   start1<-c(100,1,.1)
   start2<-c(75,2.5,0.07)
   start0<-rep(1,3)

}
