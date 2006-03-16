#  check hessian
if(!require("numDeriv"))stop("this test requires numDeriv.")

func1 <- function(x) sum(sin(x))

x <- (0:2)*2*pi/2
#x <- (0:10)*2*pi/10
#anal.h <- NEED ANALYTIC

calc.h <- hessian(func1, x) 

func2 <- function(x) cos(x)
calc.j <- jacobian(func2, x) # NOT

#if( max(abs(calc.h - calc.j)) > 1e-9) stop("hessian test FAILED")

#if( max(abs(calc.h - anal.h)) > 1e-9) stop("hessian test FAILED")

func3 <- function(x) sum(sin(x) * cos(x))

x <- (0:10)*2*pi/10
anal.h <- sin(cos(x))  # NOT 
calc.h  <- hessian(func3, x) 

#if( max(abs(calc.h - anal.h)) > 1e-9) stop("hessian test FAILED")

