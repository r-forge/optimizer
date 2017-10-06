## @knitr ##Lanczos1.prb
probname <- "##Lanczos1"
probdesc <- "

NIST/ITL StRD
Dataset Name:  Lanczos1          (Lanczos1.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 46)
               Certified Values  (lines 41 to 51)
               Data              (lines 61 to 84)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are taken from an example discussed in
               Lanczos (1956).  The data were generated to 14-digits
               of accuracy using
               f(x) = 0.0951*exp(-x) + 0.8607*exp(-3*x) 
                                     + 1.5576*exp(-5*x).


Reference:     Lanczos, C. (1956).
               Applied Analysis.
               Englewood Cliffs, NJ:  Prentice Hall, pp. 272-280.

Data:          1 Response  (y)
               1 Predictor (x)
               24 Observations
               Average Level of Difficulty
               Generated Data

Model:         Exponential Class
               6 Parameters (b1 to b6)

               y = b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x)  +  e

          Starting values                  Certified Values
 
        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   1.2         0.5           9.5100000027E-02  5.3347304234E-11
  b2 =   0.3         0.7           1.0000000001E+00  2.7473038179E-10
  b3 =   5.6         3.6           8.6070000013E-01  1.3576062225E-10
  b4 =   5.5         4.2           3.0000000002E+00  3.3308253069E-10
  b5 =   6.5         4             1.5575999998E+00  1.8815731448E-10
  b6 =   7.6         6.3           5.0000000001E+00  1.1057500538E-10

Residual Sum of Squares:                    1.4307867721E-25
Residual Standard Deviation:                8.9156129349E-14
Degrees of Freedom:                                18
Number of Observations:                            24


Data:   y                   x
       2.513400000000E+00  0.000000000000E+00
       2.044333373291E+00  5.000000000000E-02
       1.668404436564E+00  1.000000000000E-01
       1.366418021208E+00  1.500000000000E-01
       1.123232487372E+00  2.000000000000E-01
       9.268897180037E-01  2.500000000000E-01
       7.679338563728E-01  3.000000000000E-01
       6.388775523106E-01  3.500000000000E-01
       5.337835317402E-01  4.000000000000E-01
       4.479363617347E-01  4.500000000000E-01
       3.775847884350E-01  5.000000000000E-01
       3.197393199326E-01  5.500000000000E-01
       2.720130773746E-01  6.000000000000E-01
       2.324965529032E-01  6.500000000000E-01
       1.996589546065E-01  7.000000000000E-01
       1.722704126914E-01  7.500000000000E-01
       1.493405660168E-01  8.000000000000E-01
       1.300700206922E-01  8.500000000000E-01
       1.138119324644E-01  9.000000000000E-01
       1.000415587559E-01  9.500000000000E-01
       8.833209084540E-02  1.000000000000E+00
       7.833544019350E-02  1.050000000000E+00
       6.976693743449E-02  1.100000000000E+00
       6.239312536719E-02  1.150000000000E+00

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
lanczos.formula <- ( y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x)) # Lanczos


lanczos.f <- function(x) {
   res<-lanczos.res(x)
   f<-sum(res*res)
}

lanczos.res <- function(b) {
   xx<-lanczos$x # note case!
   yy<-lanczos$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   res<-b1*exp(-b2*xx) + b3*exp(-b4*xx) + b5*exp(-b6*xx) - yy
   return(res)
}

# lanczos - Jacobian
lanczos.jac <- function(b) {
stop("not defined")
   xx<-lanczos$x
   yy<-lanczos$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

lanczos.h <- function(x) {
stop("not defined")
   JJ<-lanczos.jac(x)
   H <- t(JJ) %*% JJ
   res<-lanczos.res(x)
stop("not defined")

}

lanczos.g<-function(x) {
#   stop("not defined")
   JJ<-lanczos.jac(x)
   res<-lanczos.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

lanczos.fgh<-function(x) {
   f<-lanczos.f(x)
   g<-lanczos.g(x)
   H<-lanczos.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

probdata<-function(dataname, pkgname) {
  library(pkgname, character.only=TRUE) # get parent collection
  eval(parse(text=data(list=dataname))) # and load up the data into x and y
}

lanczos.setup <- function(idxdata=NULL, idxstart=NULL) {
  if(is.null(idxdata)) {
    cat("There are three Lanczos data sets, indexed 1:3\n")     
  } else {
    if (idxdata %in% 1:3) {
      dsetname <- paste("Lanczos",idxdata,sep='')
      mypdata <- probdata(dsetname, "NISTnls")
    } else { stop("ERROR-Lanczos data set index out of range = ",idxdata)}
  }
  if(is.null(idxstart)) {
    cat("There are two Lanczos starting points, indexed 1:2, or 0 to enter from keyboard\n")     
  } else {
    if (idxstart %in% 0:2) {
      if (idxstart == 0) {
        cat("Enter the 6 elements of the starting vector:")
        for (jj in 1:6) {
          start[jj] <- readline("? ")
        }
      } else { 
        if (idxstart == 1) { 
          start<-c(1.2,0.3,5.6,5.5,6.5,7.6)
        } else { 
          start<-c(0.5,0.7,3.6,4.2,4,6.3) 
        }
        names(start)<-c("b1","b2","b3","b4","b5","b6")
      }
    } else { stop("ERROR-Lanczos starting point index out of range = ",idxdata)}
  }
  Lz <- list(mypdata=mypdata, start=start)
}



library("NISTnls", character.only=TRUE) # get parent collection
mypdata <- eval(parse(text=data(list="Lanczos1"))) # and load up the data into x and y
print(mypdata)
start1<-c(1.2,0.3,5.6,5.5,6.5,7.6)
start2<-c(0.5,0.7,3.6,4.2,4,6.3)
names(start1)<-c("b1","b2","b3","b4","b5","b6")
names(start2)<-c("b1","b2","b3","b4","b5","b6")

#lset <- lanczos.setup(1,1)
#print(lset)

library(nlsr)

st12 <- c(1.2, 0.3, 0,0,0,0)
names(st12) <- names(start1)
lanczos1nlxb12 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata, masked=c("b3","b4","b5","b6"))
print(lanczos1nlxb12)

st14 <- c(1.2,0.3,5.6,5.5,0,0)
names(st14) <- names(start1)
lanczos1nlxb14 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata, masked=c("b5","b6"))
print(lanczos1nlxb14)

st16 <- c(1.2,0.3,5.6,5.5,6.5,7.6)
names(st16) <- names(start1)
lanczos1nlxb16 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata)
print(lanczos1nlxb16)

l1start1nlxb <-  nlxb(lanczos.formula, start1, trace=TRUE, data=mypdata)
print(l1start1nlxb)
l1start2nlxb <-  nlxb(lanczos.formula, start2, trace=TRUE, data=mypdata)
print(l1start2nlxb)

