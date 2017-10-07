
## @knitr ##Candlestick
# This is file ##Candlestick.prb
probname <- "##Candlestick"
probdesc <- "
Author: John C. Nash
This problem was devised to provide a relatively straightforward
function that has an obvious minimum but that minimum is not
unique. In this function, the squared distance from the origin, r2,
is calculated. The function value is then the sum of r2 and alpha/r2,
where alpha defaults to 0.1.
To avoid divide by zero, r2 is replaced by smallval (default 1e-4)
if it is smaller than smallval.
Note that the KKT conditions are difficult to computer because the
Hessian is singular at the solution. However, most solvers have
no difficulty in finding an admissible minimum, though this will 
not be unique, as the minimum is on the hypersphere of radius
sqrt(r2) for r2 the minimum of r2 + alpha/r2.
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
##Candlestick.formula <- ( y ~ b1*x**b2 )

#- setup
# ##Candlestick.setup<-function() {
##Candlestick.df <- data.frame(x,y)

##Candlestick.start <- function(indx) {
##  ...
##  start
##}
#- Problem types will be determined from available functions / formulas
#- ?? We could have problems with quotation marks WITHIN the code.
#- lower = NA # if bounds are present, then we use bounded methods
#- upper = NA # 

#- function
##Candlestick.f <- function(x) { # ?? should remove -- this is a sumsquares problem
# candlestick function
# J C Nash 2011-2-3
cstick.f<-function(x,alpha=1e-1, smallval=1e-4){
  x<-as.vector(x)
  r2<-max(crossprod(x), smallval)
  f<-as.double(r2+alpha/r2)
  f
}

#- gradient
##Candlestick.g<-function(x) {
#-   stop('not defined')
return(gg)
}

#- hessian
##Candlestick.h <- function(x) {
cstick.g<-function(x,alpha=1e-1, smallval=1e-4){
  x<-as.vector(x)
  r2<-max(crossprod(x), smallval)
  g1<-2*x
  g2 <- (-alpha)*2*x/(r2*r2)
  g<-as.double(g1+g2)
  g
}

#- residual
##Candlestick.res <- function(b) {
#return(res)
#}

#- Jacobian
# DanielWood - Jacobian
##Candlestick.jac <- function(b) {
#return(J)
#}

#- nleq -- ?? section if appropriate


#- test
## put example calls of the function, possibly including calls to 
# optimizations and nonlinear least squares etc.



#- End Candlestick.prb   
