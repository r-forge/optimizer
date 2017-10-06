## @knitr ##Lanczos3.prb
probname <- "##Lanczos3"
probdesc <- "
NIST/ITL StRD
Dataset Name:  Lanczos3          (Lanczos3.dat)

File Format:   ASCII
Starting Values   (lines 41 to 46)
Certified Values  (lines 41 to 51)
Data              (lines 61 to 84)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are taken from an example discussed in
Lanczos (1956).  The data were generated to 5-digits
of accuracy using
f(x) = 0.0951*exp(-x) + 0.8607*exp(-3*x) 
+ 1.5576*exp(-5*x).


Reference:     Lanczos, C. (1956).
Applied Analysis.
Englewood Cliffs, NJ:  Prentice Hall, pp. 272-280.




Data:          1 Response  (y)
1 Predictor (x)
24 Observations
Lower Level of Difficulty
Generated Data

Model:         Exponential Class
6 Parameters (b1 to b6)

y = b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x)  +  e



Starting values                  Certified Values

Start 1     Start 2           Parameter     Standard Deviation
b1 =   1.2         0.5           8.6816414977E-02  1.7197908859E-02
b2 =   0.3         0.7           9.5498101505E-01  9.7041624475E-02
b3 =   5.6         3.6           8.4400777463E-01  4.1488663282E-02
b4 =   5.5         4.2           2.9515951832E+00  1.0766312506E-01
b5 =   6.5         4             1.5825685901E+00  5.8371576281E-02
b6 =   7.6         6.3           4.9863565084E+00  3.4436403035E-02

Residual Sum of Squares:                    1.6117193594E-08
Residual Standard Deviation:                2.9923229172E-05
Degrees of Freedom:                                18
Number of Observations:                            24








Data:   y           x
2.5134E+00  0.00000E+00
2.0443E+00  5.00000E-02
1.6684E+00  1.00000E-01
1.3664E+00  1.50000E-01
1.1232E+00  2.00000E-01
0.9269E+00  2.50000E-01
0.7679E+00  3.00000E-01
0.6389E+00  3.50000E-01
0.5338E+00  4.00000E-01
0.4479E+00  4.50000E-01
0.3776E+00  5.00000E-01
0.3197E+00  5.50000E-01
0.2720E+00  6.00000E-01
0.2325E+00  6.50000E-01
0.1997E+00  7.00000E-01
0.1723E+00  7.50000E-01
0.1493E+00  8.00000E-01
0.1301E+00  8.50000E-01
0.1138E+00  9.00000E-01
0.1000E+00  9.50000E-01
0.0883E+00  1.00000E+00
0.0783E+00  1.05000E+00
0.0698E+00  1.10000E+00
0.0624E+00  1.15000E+00


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
mypdata <- eval(parse(text=data(list="Lanczos3"))) # and load up the data into x and y
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
Lanczos3nlxb12 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata, masked=c("b3","b4","b5","b6"))
print(Lanczos3nlxb12)

st14 <- c(1.2,0.3,5.6,5.5,0,0)
names(st14) <- names(start1)
Lanczos3nlxb14 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata, masked=c("b5","b6"))
print(Lanczos3nlxb14)

st16 <- c(1.2,0.3,5.6,5.5,6.5,7.6)
names(st16) <- names(start1)
Lanczos3nlxb16 <- nlxb(lanczos.formula, st12, trace=TRUE, data=mypdata)
print(Lanczos3nlxb16)

l1start1nlxb <-  nlxb(lanczos.formula, start1, trace=TRUE, data=mypdata)
print(l1start1nlxb)
l1start2nlxb <-  nlxb(lanczos.formula, start2, trace=TRUE, data=mypdata)
print(l1start2nlxb)

