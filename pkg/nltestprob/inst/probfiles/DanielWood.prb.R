## @knitr DanielWood.prb
# This is file DanielWood.prb
probname <- "DanielWood"
probdesc <- "These data and model are described in Daniel and Wood (1980), and
originally published in E.S.Keeping, 'Introduction to Statistical Inference,'
Van Nostrand Company, Princeton, NJ, 1962, p. 354. The response variable is 
energy radieted from a carbon filament lamp per cm**2 per second, and the 
predictor variable is the absolute temperature of the filament in 1000 degrees Kelvin. 
NIST/ITL StRD
Dataset Name:  DanWood           (DanWood.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 42)
               Certified Values  (lines 41 to 47)
               Data              (lines 61 to 66)

Procedure:     Nonlinear Least Squares Regression

Description:   These data and model are described in Daniel and Wood
               (1980), and originally published in E.S.Keeping, 
               "Introduction to Statistical Inference," Van Nostrand
               Company, Princeton, NJ, 1962, p. 354.  The response
               variable is energy radieted from a carbon filament
               lamp per cm**2 per second, and the predictor variable
               is the absolute temperature of the filament in 1000
               degrees Kelvin.

Reference:     Daniel, C. and F. S. Wood (1980).
               Fitting Equations to Data, Second Edition. 
               New York, NY:  John Wiley and Sons, pp. 428-431.


Data:          1 Response Variable  (y = energy)
               1 Predictor Variable (x = temperature)
               6 Observations
               Lower Level of Difficulty
               Observed Data

Model:         Miscellaneous Class
               2 Parameters (b1 and b2)

               y  = b1*x**b2  +  e


 
          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   1           0.7           7.6886226176E-01  1.8281973860E-02
  b2 =   5           4             3.8604055871E+00  5.1726610913E-02
 
Residual Sum of Squares:                    4.3173084083E-03
Residual Standard Deviation:                3.2853114039E-02
Degrees of Freedom:                                4
Number of Observations:                            6 
 
Data:  y              x
      2.138E0        1.309E0
      3.421E0        1.471E0
      3.597E0        1.490E0
      4.340E0        1.565E0
      4.882E0        1.611E0
      5.660E0        1.680E0


"

#- Note: environment / list "pe" must already exist for runoptprob
## ?? maybe use pe for Problem Environment and put in pe and 
## data and other stuff.

if (exists("pe")) { 
      rm("pe")  
  }

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

  
## {stop("Environment 'pe' must exist before running problems.")}

DanielWood.formula <- ( y ~ b1*x**b2 )

#- setup
# DanielWood.setup<-function() {
y <- c( 2.138, 3.421, 3.597, 4.340, 4.882, 5.660)
x <- c( 1.309, 1.471, 1.490, 1.565, 1.611, 1.680)
DanielWood.df <- data.frame(x,y)
pe$df <- DanielWood.df
rm(x)
rm(y)
# }
print(DanielWood.df)
cat("and in pe:")
print(pe$df)

DanielWood.start <- function(indx) {
  #- indx is character string to allow for more general forms e.g., XRosenbrock
  ii <- as.numeric(indx)
  start <- NA
  if (ii == 1) {
     start = c(b1= 1, b2 = 5)
     attr(start,"fval") <- 149.7192
  }
  if (ii == 2) {
     start = c(b1 = 0.7,b2 = 4)
     attr(start,"fval") <- 0.1037647
  }
  start
}
#- Problem types will be determined from available functions / formulas
#- ?? We could have problems with quotation marks WITHIN the code.
#- lower = NA # if bounds are present, then we use bounded methods
#- upper = NA # 

#- function
DanielWood.f <- function(x) { # ?? should remove -- this is a sumsquares problem
res<-DanielWood.res(x)
pe$kf <- pe$kf + 1
f<-sum(res*res)
}

#- gradient
DanielWood.g<-function(x) {
#-   stop('not defined')
JJ<-DanielWood.jac(x)
res<-DanielWood.res(x)
gg<-as.vector(2.0*t(JJ) %*% res)
pe$kg <- pe$kg + 1
return(gg)
}

#- hessian
DanielWood.h <- function(x, DanielWood.df) {
#- THIS IS NOT COMPLETE??
  stop("DanielWood.h IS NOT COMPLETE")
  res<-DanielWood.res(x)
  JJ<-DanielWood.jac(x)
  pe$khess <- pe$khess + 1
  H <- t(JJ) %*% JJ
}


#- residual
DanielWood.res <- function(b, DanielWood.df) {
xx<-pe$df$x # case !!
yy<-pe$df$y
res <- rep(NA, length(xx))
b1<-b[1]
b2<-b[2]
res<-b1*(xx**b2) - yy
pe$kres <- pe$kres + 1
return(res)
}

#- Jacobian
# DanielWood - Jacobian
DanielWood.jac <- function(b) {
xx<-pe$df$x
yy<-pe$df$y
n<-length(b)
m<-length(xx)
b1<-b[1]
b2<-b[2]
J<-matrix(0,m,n) # define the size of the Jacobian
expr1 <- xx^b2
J[, 1] <- expr1
J[, 2] <- b1 * (expr1 * log(xx))
pe$kjac <- pe$kjac + 1
return(J)
}

# DanielWood - test
DanielWood.test <- function() {
# tests and examples of calls to DanielWood problem
y <- c( 2.138, 3.421, 3.597, 4.340, 4.882, 5.660)
x <- c( 1.309, 1.471, 1.490, 1.565, 1.611, 1.680)
DanielWood.df <- data.frame(x,y)
rm(x)
rm(y)
st1 <- DanielWood.start(1)
nls.sol <- nls(formula<-DanielWood.formula, data=DanielWood.df,
    start=st1, trace=TRUE)
print(nls.sol)
# summary(nls.sol)

# optimr
require(optimr)
cat("DanielWood.df:\n")
print(DanielWood.df)
# Note problem of passing data frame down through the functions. Need problem environment

cat("opm on DanielWood.f, start 1:\n")
print(st1)
opmdw1 <- opm(st1, DanielWood.f, gr="grcentral", method="ALL")
print(summary(opmdw1, order=value))

cat("opm on DanielWood.f, start 2:\n")
st2 <- DanielWood.start(2)
print(st2)
opmdw2 <- opm(st2, DanielWood.f, gr="grcentral", method="ALL")
print(summary(opmdw2, order=value))


}

DanielWood.test()



#- End DanielWood.prb   
