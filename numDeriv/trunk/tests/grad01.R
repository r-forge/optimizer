if(!require("numDeriv"))stop("this test requires numDeriv.")

###################################################################
#  sin test. scalar valued function with scalar arg
###################################################################

g.anal  <-  cos(pi)

g.calcR <-  grad(sin, pi, method="Richardson")

cat("max. abs. diff. between analytic and Richardson grad:")
print( max(abs(g.calcR - g.anal)), digits=16)

if(max(abs(g.calcR - g.anal)) > 1e-12)
    stop("Richardson grad scalar test FAILED")



g.calcS <-   grad(sin, pi, method="simple")

cat("max. abs. diff. between analytic and simple grad:")
print( max(abs(g.calcS - g.anal)), digits=16)

if(max(abs(g.calcS - g.anal)) > 1e-4)
     stop("simple grad scalar test FAILED")


###################################################################
#  sin test. vector argument, scalar result
###################################################################

func2a <- function(x) sum(sin(x))

x <- (0:10)*2*pi/10
g.anal  <- cos(x)

g.calcR <- grad(func2a, x, method="Richardson")

cat("max. abs. diff. between analytic and Richardson grad:")
print( max(abs(g.calcR - g.anal)), digits=16)

if(max(abs(g.calcR - g.anal)) > 1e-10)
   stop("Richardson grad scalar function vector argument test FAILED")


g.calcS <- grad(func2a, x, method="simple")

cat("max. abs. diff. between analytic and simple grad:")
print( max(abs(g.calcS - g.anal)), digits=16)

if(max(abs(g.calcS - g.anal)) > 1e-4)
    stop("simple grad scalar function vector argument test FAILED")


###################################################################
#  sin test. vector argument, vector result
###################################################################

x <- (0:10)*2*pi/10
g.anal <-  cos(x)

# THIS DOES NOT WORK FOR A FUNCTION WITH A VECTOR RESULT

g.calcR <-  grad(sin, x, method="Richardson")

cat("max. abs. diff. between analtic and Richardson grad:")
print( max(abs(g.calcR - g.anal)), digits=16)
if(max(abs(g.calcR - g.anal)) > 1e-12)
   stop("Richardson grad vector argument, vector result test FAILED")

g.calcS <-   grad(sin, x, method="simple")

cat("max. abs. diff. between analtic and simple grad:")
print( max(abs(g.calcS - g.anal)), digits=16)

if(max(abs(g.calcS - g.anal)) > 1e-4)
    stop("simple grad vector argument, vector result test FAILED")


