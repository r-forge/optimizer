## @knitr ##Ratkowsky3.prb
# This is file ##Ratkowsky3.prb.R
probname <- "##Ratkowsky3"
probdesc <- "Put your description in double quotes.
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
ratkowsky3.formula <- ( y ~ b1*x**b2 ) ??
# Optimization test function ratkowsky3
# ratkowsky3 from NISTnls
# ??ref...


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

ratkowsky3.setup<-function() {
   library(NISTnls) # get parent collection
   data(Ratkowsky3) # and load up the data into x and y
}

ratkowsky3.test<-function() {
   start1<-c(100,10,1,1)
   start2<-c(700,5,0.75,1.3)
   start0<-rep(1,4)
          
NIST<-list()
NIST$value<-8.7864049080E+03
NIST$par<-c(6.9964151270E+02,5.2771253025E+00,7.5962938329E-01,1.2792483859E+00)
NIST$ses<-c(1.6302297817E+01, 2.0828735829E+00, 1.9566123451E-01,6.8761936385E-01)

}
