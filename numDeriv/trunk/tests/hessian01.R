#  check hessian
if(!require("numDeriv"))stop("this test requires numDeriv.")


####################################################################

#   exp tests

####################################################################


x <- 0.25 * pi
print(calc.h <- hessian(sin, x) )
print(anal.h <-  sin(x+pi))
if( max(abs(calc.h - anal.h)) > 1e-8) stop("hessian test 1 FAILED")


func1 <- function(x) sum(sin(x))

x <- (0:2)*2*pi/2
#x <- (0:10)*2*pi/10
print(anal.h <-  0)
print(calc.h <- hessian(func1, x) )
if( max(abs(anal.h - calc.h)) > 1e-9) stop("hessian test FAILED")
#SOMETHING DOES NOT SEEM RIGHT NEAR BOTTOM DIAG

funcD1 <- function(x) grad(sin,x)
print(calc.j <- jacobian(funcD1, x) )
if( max(abs(calc.h - calc.j)) > 1e-4) stop("hessian test FAILED")

#if( max(abs(calc.h - anal.h)) > 1e-9) stop("hessian test FAILED")

func3 <- function(x) sum(sin(x) * cos(x))


####################################################################

#   exp tests

####################################################################

fun1e <- function(x) exp(2*x)
funD1e <- function(x) 2*exp(2*x)

x <- 1
print(anal.h <- 4*exp(2*x) )
print(calc.h  <- hessian(fun1e, x) )
if( max(abs(calc.h - anal.h)) > 1e-7) stop("hessian test FAILED")

print(calc.j <- jacobian(funD1e, x) )
if( max(abs(calc.j - anal.h)) > 1e-7) stop("hessian test FAILED")

fun1e <- function(x) sum(exp(2*x))
funD1e <- function(x) 2*exp(2*x)
x <- c(1,3,5)
print(anal.h <- diag(4*exp(2*x)) )
print(calc.h  <- hessian(fun1e, x) )
if( max(abs(calc.h - anal.h)) > 1e-5) stop("hessian test FAILED")

print(calc.j <- jacobian(funD1e, x) )
if( max(abs(calc.j - anal.h)) > 1e-6) stop("hessian test FAILED")
