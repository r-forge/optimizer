## @knitr ##FletcherPowellHelicalValley.prb
# This is file ##FletcherPowellHelicalValley.prb
probname <- "##FletcherPowellHelicalValley"
probdesc <- "From Fletcher and Powell 1963:
  A function with a steep sided helical valley. This function
     f(x1, x2, x3) = 100{[x3 - 10 theta(x1, x2)]^2 
                     + [r(x1, x2) - 1]^2 } + x3^2

where

2 pi theta(x1, x2) = arctan(x2/x1) for x1 > 0
                   = pi + arctan (x2/x1) for x1 < 0

and

   r(x1, x2) = (x1^2 + x2^2)^(0.5)

has a helical valley in the x3 direction with pitch 10 and
radius 1. It is only considered for

      — pi/2 < 2 pi theta < 3 pi/2

that is

           —2.5 < x3 < 7.5

JN: ?? Thy does this follow. theta is about x1 and x2, not x3.
What does pitch 10 and radius 1 mean?
"


#*****************************************************************************
## fphelical evaluates the Fletcher-Powell helical valley function.

# JN 170920 -- note that Burkardt helical.m has errors!!!
#
#  Parameters:
#
#    Input, real x of length 3, the argument.
#
#    Output, real f, the value of the function.
#
# atan2(y, x) = atan(y/x)

## jbfphelical.f <- function(x) { # from Burkardt helical.m
##    if ( x[1] > 0.0 ) {
##       theta <- atan2 ( x[2], x[1] ) / 2.0 / pi
##    } else if ( x[1] < 0.0 ) {
##       theta <- 0.5 + atan2 ( x[2], x[1] ) / 2.0 / pi
##    } else if ( x[1] == 0.0 ) { theta <- 0.25 }
##           else stop("fphelical.f -- impossible condition!")
##   fx1 = x[3] - 10.0 * theta
##   fx2 = sqrt ( x[1] * x[1] + x[2] * x[2] )
##   fx3 = x[3]
##   fx = 100.0 * fx1 * fx1 + fx2 * fx2 +  fx3 * fx3;
###  Note: this is a sum of squares -- should be able to produce nls form
##   fx
## }


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
##FletcherPowellHelicalValley.formula <- ( y ~ b1*x**b2 )

#- setup
# ##FletcherPowellHelicalValley.setup<-function() {
y <- c( 2.138, 3.421, 3.597, 4.340, 4.882, 5.660) # example only
x <- c( 1.309, 1.471, 1.490, 1.565, 1.611, 1.680)
##FletcherPowellHelicalValley.df <- data.frame(x,y)
rm(x)
rm(y)

##FletcherPowellHelicalValley.start <- function(indx) {
#  ...

#  start
#}
#- Problem types will be determined from available functions / formulas
#- ?? We could have problems with quotation marks WITHIN the code.
#- lower = NA # if bounds are present, then we use bounded methods
#- upper = NA # 

#- function
##FletcherPowellHelicalValley.f <- function(x) { # ?? should remove -- this is a 

fphelical.f <- function(x) {#  From Fletcher-Powell 1963 paper on Davidon method
   if (x[1] == 0) {theta <- sign(x[2])*1e+20 } # safety setting
   else { if ( x[1] > 0.0 ) { theta <- atan(x[2]/x[1])/(2*pi) }
           else { theta <- 0.5 + atan(x[2]/x[1])/(2.0*pi) }
   }
   r <- sqrt(x[1] * x[1] + x[2] * x[2])
   f <- 100 * ( (x[3] - 10*theta)^2 + (r-1)^2 ) + x[3]^2
}

#- gradient
##FletcherPowellHelicalValley.g<-function(x) {
fphelical.g <- function(x) {
   if (x[1] == 0) {theta <- sign(x[2])*1e+20 } # safety setting
   else { if ( x[1] > 0.0 ) { theta <- atan(x[2]/x[1])/(2*pi) }
           else { theta <- 0.5 + atan(x[2]/x[1])/(2.0*pi) }
   }
   r <- sqrt(x[1] * x[1] + x[2] * x[2])
  if (x[1] == 0) { gt2 <- 0 }
  else {  gt2 <- (1/x[1])/((1+(x[2]/x[1])^2)*(2*pi)) }

  if (x[1] == 0) { gt1 <- 0 }
  else {  gt1 <- (-x[2]/(x[1]^2))/((1+(x[2]/x[1])^2)*(2*pi)) }

  r1 <- x[1]/r
  r2 <- x[2]/r   
  g1 <- 100*( -20*(x[3]-10*theta)*gt1 + 2*(r-1)*r1)
  g2 <- 100*( -20*(x[3]-10*theta)*gt2 + 2*(r-1)*r2)

  g3 <- 202*x[3] - 2000*theta
  gg <- c(g1, g2, g3)
  gg
}

#- hessian
##FletcherPowellHelicalValley.h <- function(x) {
#- THIS IS NOT COMPLETE??


#- residual
##FletcherPowellHelicalValley.res <- function(b) {

#- Jacobian
##FletcherPowellHelicalValley.jac <- function(b) {

#- nleq -- ?? section if appropriate


#- test
## put example calls of the function, possibly including calls to 
# optimizations and nonlinear least squares etc.

library(numDeriv)

fphelical.test <- function() {
  x0 <- c(-1, 0, 0)
  xstar <- c(1, 0, 0)
  
  f0fp <- fphelical.f(x0)
  cat("x0=")
  print(x0)
  cat("fphelical(x0)=")
  print(f0fp)
  cat("\n")
  ga0 <- fphelical.g(x0)
  gn0 <- grad(fphelical.f, x0)
  cat("gradient at x0:")
  print(ga0)
  cat("numerical  gn0:")
  print(gn0)
  cat("max(abs(diff))=", max(abs(ga0-gn0)),"\n")
  
  cat("xstar=")
  print(xstar)
  cat("fphelical(xstar)=")  
  fstarfp <- fphelical.f(xstar)
  print(fstarfp)
  gas <- fphelical.g(xstar)
  gns <- grad(fphelical.f, xstar)
  cat("gradient at xstar:")
  print(gas)
  cat("numerical  gnstar:")
  print(gns)
  cat("max(abs(diff))=", max(abs(gas-gns)),"\n")
}

fphelical.test() # To verify that function is working correctly

## Examples
library(optimr)
x0 <- c(-1, 0, 0)
fphall <- opm(x0, fphelical.f, fphelical.g, method="ALL")

print(summary(fphall, order=value))



#- End FletcherPowellHelicalValley.prb   



