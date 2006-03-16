#  check jacobian

if(!require("numDeriv"))stop("this test requires numDeriv.")


func2 <- function(x) c(sin(x), cos(x))

x <- (1:2)*2*pi/2
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test FAILED")

x <- (0:1)*2*pi
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test FAILED")


x <- (0:10)*2*pi/10
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test FAILED")
