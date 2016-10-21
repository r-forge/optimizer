## ------------------------------------------------------------------------

dx2x <- deriv(~ x^2, "x") 
dx2x
mode(dx2x)
str(dx2x)
x <- -1:2
eval(dx2x) # This is evaluated at -1, 0, 1, 2, with the result in the "gradient" attribute
# Note that we cannot (easily) differentiate this again.
firstd <- attr(dx2x,"gradient")
str(firstd)
d2x2x <- try(deriv(firstd, "x"))
str(d2x2x)
### Now try D()
Dx2x <- D(expression(x^2), "x")
Dx2x
x <- -1:2
eval(Dx2x)
# We can differentiate aggain
D2x2x <- D(Dx2x,"x")
D2x2x
eval(D2x2x) # But we don't get a vector -- could be an issue in gradients/Jacobians
# Note how we handle an expression stored in a string via parse(text=  ))
sx2 <- "x^2"
sDx2x <- D(parse(text=sx2), "x")
sDx2x
# But watch out! The following "seems" to work, but the answer is not as intended
Dx2xx <- D((x^2), "x")
Dx2xx
eval(Dx2xx)

## Something 'tougher':
trig.exp <- expression(sin(cos(x + y^2)))
( D.sc <- D(trig.exp, "x") )
all.equal(D(trig.exp[[1]], "x"), D.sc)

( dxy <- deriv(trig.exp, c("x", "y")) )
y <- 1
eval(dxy)
eval(D.sc)

## function returned:
deriv((y ~ sin(cos(x) * y)), c("x","y"), func = TRUE)
##??## Surely there is an error, since docs say no lhs!

## function with defaulted arguments:
(fx <- deriv(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),
             function(b0, b1, th, x = 1:7){} ) )
fx(2, 3, 4)

## First derivative

D(expression(x^2), "x")
stopifnot(D(as.name("x"), "x") == 1)

## Higher derivatives
myd3 <- deriv3(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),
     c("b0", "b1", "th", "x") )
myd3(2,3,4, x=1:7)


## Higher derivatives:
DD <- function(expr, name, order = 1) {
   if(order < 1) stop("'order' must be >= 1")
   if(order == 1) D(expr, name)
   else DD(D(expr, name), name, order - 1)
}
DD(expression(sin(x^2)), "x", 3)
## showing the limits of the internal "simplify()" :
## -sin(x^2) * (2 * x) * 2 + ((cos(x^2) * (2 * x) * (2 * x) + sin(x^2) *
##    2) * (2 * x) + sin(x^2) * (2 * x) * 2)


## ------------------------------------------------------------------------
require(nlsr)

# Examples

newDeriv()
newDeriv(sin(x))
Deriv(sin(x+y), "x")

f <- function(x) x^2 # a "new" function
newDeriv(f(x), 2*x*D(x)) # and its derivative definition
Deriv(f(abs(x)), "x") # and the derivative of the function of abs(x)

# try fnDeriv
myfd <- fnDeriv(f(abs(x)), "x") # and the derivative of the function of abs(x)
myfd

# a different example
Deriv(pnorm(x, sd=2, log = TRUE), "x")











# require(stats)
# require(Deriv)
# require(Ryacas)

# Various derivatives 

new <- fnDeriv(quote(1 + x + y), c("x", "y"))
old <- deriv(quote(1 + x + y), c("x", "y"))
print(new)
# Following generates a very long line on output of knitr (for markdown)
class(new)
str(new)
as.expression(new)
print(old)
class(old)
str(old)

## ------------------------------------------------------------------------
fnfromnew <- function(x,y){
    .value <- 1 + x + y
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("x", 
    "y")))
    .grad[, "x"] <- 1
    .grad[, "y"] <- 1
    attr(.value, "gradient") <- .grad
    .value
}

print(fnfromnew(3,5))


## ------------------------------------------------------------------------
x <- NA
y <- Inf
print(eval(new))
print(eval(old))

## ------------------------------------------------------------------------
safeD <- function(obj, var) {
   # safeguarded D() function for symbolic derivs
   if (! is.character(var) ) stop("The variable var MUST be character type")
   if (is.character(obj) ) {
       eobj <- parse(text=obj)
       result <- D(eobj, var)
   } else {
       result <- D(obj, var)
   }
}

lxy2 <- expression(log(x+y^2))
clxy2 <- "log(x+y^2)"
try(print(D(clxy2, "y")))
print(try(D(lxy2, "y")))
print(safeD(clxy2, "y"))
print(safeD(lxy2, "y"))

## ------------------------------------------------------------------------
## Try different ways to supply the log function
aDeriv <- Deriv(log(x), "x")
class(aDeriv)
aDeriv
aderiv <- try(deriv( ~ log(x), "x"))
class(aderiv)
aderiv
aD <- D(expression(log(x)), "x")
class(aD)
aD
cat("but \n")
try(D( "~ log(x)", "x")) # fails -- gives NA rather than expected answer due to quotes
try(D( ~ log(x), "x"))
interm <- ~ log(x)
interm
class(interm)
interme <- as.expression(interm)
class(interme)
try(D(interme, "x"))
try(deriv(interme, "x"))
try(deriv(interm, "x"))


Deriv(log(x, base=3), "x" ) # OK
try(D(expression(log(x, base=3)), "x" )) # fails - only single-argument calls supported
try(deriv(~ log(x, base=3), "x" )) # fails - only single-argument calls supported
try(deriv(expression(log(x, base=3)), "x" )) # fails - only single-argument calls supported
try(deriv3(expression(log(x, base=3)), "x" )) # fails - only single-argument calls supported
fnDeriv(quote(log(x, base=3)), "x" )

Deriv(exp(x), "x")
D(expression(exp(x)), "x") # OK
deriv(~exp(x), "x") # OK, but much more complicated
fnDeriv(quote(exp(x)), "x")

Deriv(sin(x), "x")
D(expression(sin(x)), "x")
deriv(~sin(x), "x")
fnDeriv(quote(sin(x)), "x")

Deriv(cos(x), "x")
D(expression(cos(x)), "x")
deriv(~ cos(x), "x")
fnDeriv(quote(cos(x)), "x")

Deriv(tan(x), "x")
D(expression(tan(x)), "x")
deriv(~ tan(x), "x")
fnDeriv(quote(tan(x)), "x")

Deriv(sinh(x), "x")
D(expression(sinh(x)), "x")
deriv(~sinh(x), "x")
fnDeriv(quote(sinh(x)), "x")

Deriv(cosh(x), "x")
D(expression(cosh(x)), "x")
deriv(~cosh(x), "x")
fnDeriv(quote(cosh(x)), "x")

Deriv(sqrt(x), "x")
D(expression(sqrt(x)), "x")
deriv(~sqrt(x), "x")
fnDeriv(quote(sqrt(x)), "x")

Deriv(pnorm(q), "q")
D(expression(pnorm(q)), "q")
deriv(~pnorm(q), "q")
fnDeriv(quote(pnorm(q)), "q")

Deriv(dnorm(x, mean), "mean")
D(expression(dnorm(x, mean)), "mean")
deriv(~dnorm(x, mean), "mean")
fnDeriv(quote(dnorm(x, mean)), "mean")

Deriv(asin(x), "x")
D(expression(asin(x)), "x")
deriv(~asin(x), "x")
fnDeriv(quote(asin(x)), "x")

Deriv(acos(x), "x")
D(expression(acos(x)), "x")
deriv(~acos(x), "x")
fnDeriv(quote(acos(x)), "x")

Deriv(atan(x), "x")
D(expression(atan(x)), "x")
deriv(~atan(x), "x")
fnDeriv(quote(atan(x)), "x")

Deriv(gamma(x), "x")
D(expression(gamma(x)), "x")
deriv(~gamma(x), "x")
fnDeriv(quote(gamma(x)), "x")

Deriv(lgamma(x), "x")
D(expression(lgamma(x)), "x")
deriv(~lgamma(x), "x")
fnDeriv(quote(lgamma(x)), "x")

Deriv(digamma(x), "x")
D(expression(digamma(x)), "x")
deriv(~digamma(x), "x")
fnDeriv(quote(digamma(x)), "x")

Deriv(trigamma(x), "x")
D(expression(trigamma(x)), "x")
deriv(~trigamma(x), "x")
fnDeriv(quote(trigamma(x)), "x")

Deriv(psigamma(x, deriv = 5), "x")
D(expression(psigamma(x, deriv = 5)), "x")
deriv(~psigamma(x, deriv = 5), "x")
fnDeriv(quote(psigamma(x, deriv = 5)), "x")

Deriv(x*y, "x")
D(expression(x*y), "x")
deriv(~x*y, "x")
fnDeriv(quote(x*y), "x")

Deriv(x/y, "x")
D(expression(x/y), "x")
deriv(~x/y, "x")
fnDeriv(quote(x/y), "x")

Deriv(x^y, "x")
D(expression(x^y), "x")
deriv(~x^y, "x")
fnDeriv(quote(x^y), "x")

Deriv((x), "x")
D(expression((x)), "x")
deriv(~(x), "x")
fnDeriv(quote((x)), "x")

Deriv(+x, "x")
D(expression(+x), "x")
deriv(~ +x, "x")
fnDeriv(quote(+x), "x")

Deriv(-x, "x")
D(expression(- x), "x")
deriv(~ -x, "x")
fnDeriv(quote(-x), "x")

Deriv(abs(x), "x")
try(D(expression(abs(x)), "x")) # 'abs' not in derivatives table
try(deriv(~ abs(x), "x"))
fnDeriv(quote(abs(x)), "x")

Deriv(sign(x), "x")
try(D(expression(sign(x)), "x")) # 'sign' not in derivatives table
try(deriv(~ sign(x), "x"))
fnDeriv(quote(sign(x)), "x")

## ------------------------------------------------------------------------
# Various simplifications

Simplify(quote(+(a+b)))
Simplify(quote(-5))
Simplify(quote(--(a+b)))

Simplify(quote(exp(log(a+b))))
Simplify(quote(exp(1)))

Simplify(quote(log(exp(a+b))))
Simplify(quote(log(1)))

Simplify(quote(!TRUE))
Simplify(quote(!FALSE))

Simplify(quote((a+b)))

Simplify(quote(a + b + 0))
Simplify(quote(0 + a + b))
Simplify(quote((a+b) + (a+b)))
Simplify(quote(1 + 4))

Simplify(quote(a + b - 0))
Simplify(quote(0 - a - b))
Simplify(quote((a+b) - (a+b)))
Simplify(quote(5 - 3))

Simplify(quote(0*(a+b)))
Simplify(quote((a+b)*0))
Simplify(quote(1L * (a+b)))
Simplify(quote((a+b) * 1))
Simplify(quote((-1)*(a+b)))
Simplify(quote((a+b)*(-1)))
Simplify(quote(2*5))

Simplify(quote((a+b) / 1))
Simplify(quote((a+b) / (-1)))
Simplify(quote(0/(a+b)))
Simplify(quote(1/3))

Simplify(quote((a+b) ^ 1))
Simplify(quote(2^10))

Simplify(quote(log(exp(a), 3)))

Simplify(quote(FALSE && b))
Simplify(quote(a && TRUE))
Simplify(quote(TRUE && b))

Simplify(quote(a || TRUE))
Simplify(quote(FALSE || b))
Simplify(quote(a || FALSE))

Simplify(quote(if (TRUE) a+b))
Simplify(quote(if (FALSE) a+b))

Simplify(quote(if (TRUE) a+b else a*b))
Simplify(quote(if (FALSE) a+b else a*b))
Simplify(quote(if (cond) a+b else a+b))

# This one was wrong...
Simplify(quote(--(a+b)))


## ------------------------------------------------------------------------
require(nlsr)
dnlsr <- Deriv(sin(x+y), "x")
print(dnlsr)
class(dnlsr)

require(Ryacas)
x <- Sym("x")
y <- Sym("y")
dryacas <- deriv(sin(x+y), x)
print(dryacas)
class(dryacas)


## ------------------------------------------------------------------------
require(nlsr)
dlogx <- nlsr::Deriv(log(x), "x")
str(dlogx)
print(dlogx)

## ------------------------------------------------------------------------
dlogxs <- nlsr::Deriv(expression(log(x)), "x", do_substitute=FALSE)
str(dlogxs)
print(dlogxs)
cat(as.character(dlogxs), "\n")
fne <- expression(log(x))
dlogxe <- Deriv(fne, "x", do_substitute=FALSE)
str(dlogxe)
print(dlogxe)


# base R
dblogx <- D(expression(log(x)), "x")
str(dblogx)
print(dblogx)

require(Deriv)
ddlogx <- Deriv::Deriv(expression(log(x)), "x")
str(ddlogx)
print(ddlogx)
cat(as.character(ddlogx), "\n")
ddlogxf <- ~ ddlogx
str(ddlogxf)

## ------------------------------------------------------------------------
zzz <- expression(y[3]*r1 + r2)
try(deriv(zzz,c("r1","r2")))
require(nlsr)
try(nlsr::Deriv(zzz, c("r1","r2")))
try(fnDeriv(zzz, c("r1","r2")))
newDeriv(`[`(x,y), stop("no derivative when indexing"))
try(nlsr::Deriv(zzz, c("r1","r2")))
try(nlsr::fnDeriv(zzz, c("r1","r2")))

## ------------------------------------------------------------------------
try(nlsr::Deriv(zzz, "y[3]"))
try(nlsr::Deriv(y3*r1+r2,"y3"))
try(nlsr::Deriv(y[3]*r1+r2,"y[3]"))

