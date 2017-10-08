## @knitr ##ENSO.prb
# This is file ##ENSO.prb
probname <- "##ENSO"
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
##ENSO.formula <- ( y ~ b1*x**b2 )

#- setup

library(NISTnls, character.only=TRUE)
mypdata <- eval(parse(text=data("ENSO")))# Optimization test function enso
# enso from NISTnls
# ??ref...


enso.f <- function(x) {
   res<-enso.res(x)
   f<-sum(res*res)
}

enso.res <- function(b) {
# NOTE: could benefit from some sort of constraint to avoid equal parameters in trig args.
   xx<-ENSO$x # note case!
   yy<-ENSO$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   b7<-b[7]
   b8<-b[8]
   b9<-b[9]
   res<-b1 + b2*cos( 2*pi*xx/12 ) + b3*sin( 2*pi*xx/12 ) + b5*cos( 2*pi*xx/b4 ) + b6*sin( 2*pi*xx/b4 ) + b8*cos( 2*pi*xx/b7 ) + b9*sin( 2*pi*xx/b7 )  - yy
   return(res)
}

# enso - Jacobian
enso.jac <- function(b) {
stop("not defined")
   xx<-enso$x
   yy<-enso$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   J<-matrix(0,m,n) # define the size of the Jacobian
    return(J)
}

enso.h <- function(x) {
stop("not defined")
   JJ<-enso.jac(x)
   H <- t(JJ) %*% JJ
   res<-enso.res(x)
stop("not defined")

}

enso.g<-function(x) {
#   stop("not defined")
   JJ<-enso.jac(x)
   res<-enso.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

enso.fgh<-function(x) {
   f<-enso.f(x)
   g<-enso.g(x)
   H<-enso.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

enso.setup<-function() {
   library(NISTnls) # get parent collection
   data(ENSO) # and load up the data into x and y
}

enso.test<-function() {
}   
