## @knitr ##Thurber.prb
# This is file ##Thurber.prb
rm(list=ls())
probname <- "##Thurber"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Thurber           (Thurber.dat)

File Format:   ASCII
Starting Values   (lines 41 to 47)
Certified Values  (lines 41 to 52)
Data              (lines 61 to 97)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study involving
semiconductor electron mobility.  The response 
variable is a measure of electron mobility, and the 
predictor variable is the natural log of the density.


Reference:     Thurber, R., NIST (197?).  
Semiconductor electron mobility modeling.






Data:          1 Response Variable  (y = electron mobility)
1 Predictor Variable (x = log[density])
37 Observations
Higher Level of Difficulty
Observed Data

Model:         Rational Class (cubic/cubic)
7 Parameters (b1 to b7)

y = (b1 + b2*x + b3*x**2 + b4*x**3) / 
(1 + b5*x + b6*x**2 + b7*x**3)  +  e


Starting Values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   1000        1300          1.2881396800E+03  4.6647963344E+00
b2 =   1000        1500          1.4910792535E+03  3.9571156086E+01
b3 =    400         500          5.8323836877E+02  2.8698696102E+01
b4 =     40          75          7.5416644291E+01  5.5675370270E+00
b5 =      0.7         1          9.6629502864E-01  3.1333340687E-02
b6 =      0.3         0.4        3.9797285797E-01  1.4984928198E-02
b7 =      0.03        0.05       4.9727297349E-02  6.5842344623E-03

Residual Sum of Squares:                    5.6427082397E+03
Residual Standard Deviation:                1.3714600784E+01
Degrees of Freedom:                                30
Number of Observations:                            37







Data:   y             x
80.574E0      -3.067E0
84.248E0      -2.981E0
87.264E0      -2.921E0
87.195E0      -2.912E0
89.076E0      -2.840E0
89.608E0      -2.797E0
89.868E0      -2.702E0
90.101E0      -2.699E0
92.405E0      -2.633E0
95.854E0      -2.481E0
100.696E0      -2.363E0
101.060E0      -2.322E0
401.672E0      -1.501E0
390.724E0      -1.460E0
567.534E0      -1.274E0
635.316E0      -1.212E0
733.054E0      -1.100E0
759.087E0      -1.046E0
894.206E0      -0.915E0
990.785E0      -0.714E0
1090.109E0      -0.566E0
1080.914E0      -0.545E0
1122.643E0      -0.400E0
1178.351E0      -0.309E0
1260.531E0      -0.109E0
1273.514E0      -0.103E0
1288.339E0       0.010E0
1327.543E0       0.119E0
1353.863E0       0.377E0
1414.509E0       0.790E0
1425.208E0       0.963E0
1421.384E0       1.006E0
1442.962E0       1.115E0
1464.350E0       1.572E0
1468.705E0       1.841E0
1447.894E0       2.047E0
1457.628E0       2.200E0

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
Thurber.formula <- ( y ~ (b1 + b2*x + b3*x**2 + b4*x**3) / (1 + b5*x + b6*x**2 + b7*x**3) )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Thurber")))
# Optimization test function thurber
# thurber from NISTnls
# ??ref...


thurber.f <- function(x) {
   res<-thurber.res(x)
   f<-sum(res*res)
}

thurber.res <- function(b) {
   xx<-Thurber$x
   yy<-Thurber$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   res<-(b1 + b2*xx + b3*xx**2 + b4*xx**3) / (1 + b5*xx + b6*xx**2 + b7*xx**3) - yy
   return(res)
}

# thurber - Jacobian
thurber.jac <- function(b) {
stop("not defined")
   xx<-thurber$x
   yy<-thurber$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

thurber.h <- function(x) {
stop("not defined")
   JJ<-thurber.jac(x)
   H <- t(JJ) %*% JJ
   res<-thurber.res(x)
stop("not defined")

}

thurber.g<-function(x) {
#   stop("not defined")
   JJ<-thurber.jac(x)
   res<-thurber.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

thurber.fgh<-function(x) {
   f<-thurber.f(x)
   g<-thurber.g(x)
   H<-thurber.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}


start1<-c(1000,1000,400,40,0.7,0.3,0.03)
names(start1) <- c("b1","b2","b3","b4","b5","b6","b7")

start2<-c(1300,1500,500,75,1,0.4,0.05)
names(start2) <- c("b1","b2","b3","b4","b5","b6","b7")


start0<-rep(1,7)
names(start0) <- c("b1","b2","b3","b4","b5","b6","b7")


NISTThurber<-list()
NISTThurber$value<-5.6427082397E+03
NISTThurber$par<-c(1.2881396800E+03, 1.4910792535E+03, 5.8323836877E+02, 7.5416644291E+01, 9.6629502864E-01, 3.9797285797E-01,  4.9727297349E-02)
NISTThurber$ses<-c( 4.6647963344E+00, 3.9571156086E+01,  2.8698696102E+01,  5.5675370270E+00,  3.1333340687E-02, 1.4984928198E-02,  6.5842344623E-03)

   

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Thurber)
Try(fm1 <- nls(y ~ (b1+x*(b2+x*(b3+b4*x))) / (1+x*(b5+x*(b6+x*b7))),
               data = Thurber, trace = TRUE,
               start = c(b1 = 1000, b2 = 1000, b3 = 400, b4 = 40,
                         b5 = 0.7, b6 = 0.3, b7 = 0.03)))
Try(fm1a <- nls(y ~ (b1+x*(b2+x*(b3+b4*x))) / (1+x*(b5+x*(b6+x*b7))),
                data = Thurber, trace = TRUE, alg = "port", 
                start = c(b1 = 1000, b2 = 1000, b3 = 400, b4 = 40,
                          b5 = 0.7, b6 = 0.3, b7 = 0.03)))
Try(fm2 <- nls(y ~ (b1+x*(b2+x*(b3+b4*x))) / (1+x*(b5+x*(b6+x*b7))),
               data = Thurber, trace = TRUE,
               start = c(b1 = 1300, b2 = 1500, b3 = 500, b4 = 75,
                         b5 = 1, b6 = 0.4, b7 = 0.05)))
Try(fm2a <- nls(y ~ (b1+x*(b2+x*(b3+b4*x))) / (1+x*(b5+x*(b6+x*b7))),
                data = Thurber, trace = TRUE, alg = "port", 
                start = c(b1 = 1300, b2 = 1500, b3 = 500, b4 = 75,
                          b5 = 1, b6 = 0.4, b7 = 0.05)))
Try(fm3 <- nls(y ~ outer(x, 0:3, "^")/(1+x*(b5+x*(b6+x*b7))),
               data = Thurber, trace = TRUE,
               start = c(b5 = 0.7, b6 = 0.3, b7 = 0.03), alg = "plinear"))
Try(fm4 <- nls(y ~ outer(x, 0:3, "^")/(1+x*(b5+x*(b6+x*b7))),
               data = Thurber, trace = TRUE,
               start = c(b5 = 1, b6 = 0.4, b7 = 0.05), alg = "plinear"))

library(nlsr)

Thurbernlxb0 <- nlxb(start=start0, formula=Thurber.formula, data=mypdata, trace=TRUE)
print(Thurbernlxb0)

Thurbernlxb1 <- nlxb(start=start1, formula=Thurber.formula, data=mypdata, trace=TRUE)
print(Thurbernlxb1)

Thurbernlxb2 <- nlxb(start=start2, formula=Thurber.formula, data=mypdata, trace=TRUE)
print(Thurbernlxb2)

library(optimr)

mset  <- "ALL"
Thurberopm0fwd <- opm(start0, thurber.f, "grfwd", method=mset)
summary(Thurberopm0fwd, order=value)

Thurberopm1fwd <- opm(start1, thurber.f, "grfwd", method=mset)
summary(Thurberopm1fwd, order=value)

Thurberopm2fwd <- opm(start2, thurber.f, "grfwd", method=mset)
summary(Thurberopm2fwd, order=value)

