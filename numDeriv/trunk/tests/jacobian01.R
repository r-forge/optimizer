#  check jacobian

if(!require("numDeriv"))stop("this test requires numDeriv.")

x <- pi
j.calc <- jacobian(sin, x)
print(j.calc)
str(j.calc )
if( max(abs(j.calc - cos(x))) > 1e-10) stop("jacobian matrix test 1 FAILED")

x <- (1:2)*2*pi/2
j.calc <- jacobian(sin, x)
if( max(abs(j.calc -  diag(cos(x)))) > 1e-10) stop("jacobian matrix test 2 FAILED")

func2 <- function(x) c(sin(x), cos(x))

x <- (1:2)*2*pi/2
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test 3 FAILED")

x <- (0:1)*2*pi
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test 4 FAILED")


x <- (0:10)*2*pi/10
j.calc <- jacobian(func2, x)
if( max(abs(j.calc - 
 rbind(diag(cos(x)), diag(-sin(x))))) > 1e-10) stop("jacobian matrix test 5 FAILED")
