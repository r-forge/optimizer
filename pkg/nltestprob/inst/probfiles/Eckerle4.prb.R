# Optimization test function eckerle4
# eckerle4 from NISTnls
# ??ref...
## @knitr ##YourProblem.prb
# This is file ##YourProblem.prb
probname <- "eckerle4"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Eckerle4          (Eckerle4.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 43)
               Certified Values  (lines 41 to 48)
               Data              (lines 61 to 95)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study involving
               circular interference transmittance.  The response
               variable is transmittance, and the predictor variable
               is wavelength.


Reference:     Eckerle, K., NIST (197?).  
               Circular Interference Transmittance Study.

Data:          1 Response Variable  (y = transmittance)
               1 Predictor Variable (x = wavelength)
               35 Observations
               Higher Level of Difficulty
               Observed Data

Model:         Exponential Class
               3 Parameters (b1 to b3)

               y = (b1/b2) * exp[-0.5*((x-b3)/b2)**2]  +  e



          Starting values                  Certified Values
 
        Start 1     Start 2           Parameter     Standard Deviation
  b1 =     1           1.5         1.5543827178E+00  1.5408051163E-02
  b2 =    10           5           4.0888321754E+00  4.6803020753E-02
  b3 =   500         450           4.5154121844E+02  4.6800518816E-02

Residual Sum of Squares:                    1.4635887487E-03
Residual Standard Deviation:                6.7629245447E-03
Degrees of Freedom:                                32
Number of Observations:                            35

Data:  y                x
      0.0001575E0    400.000000E0
      0.0001699E0    405.000000E0
      0.0002350E0    410.000000E0
      0.0003102E0    415.000000E0
      0.0004917E0    420.000000E0
      0.0008710E0    425.000000E0
      0.0017418E0    430.000000E0
      0.0046400E0    435.000000E0
      0.0065895E0    436.500000E0
      0.0097302E0    438.000000E0
      0.0149002E0    439.500000E0
      0.0237310E0    441.000000E0
      0.0401683E0    442.500000E0
      0.0712559E0    444.000000E0
      0.1264458E0    445.500000E0
      0.2073413E0    447.000000E0
      0.2902366E0    448.500000E0
      0.3445623E0    450.000000E0
      0.3698049E0    451.500000E0
      0.3668534E0    453.000000E0
      0.3106727E0    454.500000E0
      0.2078154E0    456.000000E0
      0.1164354E0    457.500000E0
      0.0616764E0    459.000000E0
      0.0337200E0    460.500000E0
      0.0194023E0    462.000000E0
      0.0117831E0    463.500000E0
      0.0074357E0    465.000000E0
      0.0022732E0    470.000000E0
      0.0008800E0    475.000000E0
      0.0004579E0    480.000000E0
      0.0002345E0    485.000000E0
      0.0001586E0    490.000000E0
      0.0001143E0    495.000000E0
      0.0000710E0    500.000000E0
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

eckerle4.formula <- ( y ~ (b1/b2) * exp(-0.5*((x-b3)/b2)**2))

eckerle4.f <- function(x) {
   res<-eckerle4.res(x)
   f<-sum(res*res)
}

eckerle4.res <- function(b) {
   xx<-Eckerle4$x # note case!
   yy<-Eckerle4$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-(b1/b2) * exp(-0.5*((xx-b3)/b2)**2) - yy
   return(res)
}

# eckerle4 - Jacobian
eckerle4.jac <- function(b) {
stop("not defined")
   xx<-eckerle4$x
   yy<-eckerle4$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

eckerle4.h <- function(x) {
stop("not defined")
   JJ<-eckerle4.jac(x)
   H <- t(JJ) %*% JJ
   res<-eckerle4.res(x)
stop("not defined")

}

eckerle4.g<-function(x) {
#   stop("not defined")
   JJ<-eckerle4.jac(x)
   res<-eckerle4.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

eckerle4.fgh<-function(x) {
   f<-eckerle4.f(x)
   g<-eckerle4.g(x)
   H<-eckerle4.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

probdata<-function(dataname, pkgname) {
  library(pkgname, character.only=TRUE) # get parent collection
  eval(parse(text=data(list=dataname))) # and load up the data into x and y
}

eckerle4.setup<-function() {
   library(NISTnls) # get parent collection
#   data(Eckerle4) # and load up the data into x and y
   dsetname<-"Eckerle4"
   mypdata <- probdata(dsetname, "NISTnls")
}

eckerle4.test<-function() {
# Currently empty
}   


mypdata <- eckerle4.setup()

library(nlsr)
x1 <- c(1, 10, 500)
names(x1)<-c("b1","b2","b3")
x2 <- c(1.5, 5, 450)
names(x2)<-c("b1","b2","b3")

ecknlsr1 <- nlxb(eckerle4.formula, start=x1, data=mypdata, trace=TRUE)
# summary(ecknlsr1)
print(ecknlsr1)

ecknlsr2 <- nlxb(eckerle4.formula, start=x2, data=mypdata, trace=TRUE)
# summary(ecknlsr1)
print(ecknlsr2)
