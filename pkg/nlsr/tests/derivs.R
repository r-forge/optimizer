library(nlsr)

# Various derivatives 

fnDeriv(quote(1 + x + y), c("x", "y"))

Deriv(log(x), "x")
Deriv(log(x, base=3), "x" )

Deriv(exp(x), "x")
Deriv(sin(x), "x")
Deriv(cos(x), "x")
Deriv(tan(x), "x")
Deriv(sinh(x), "x")
Deriv(cosh(x), "x")
Deriv(sqrt(x), "x")
Deriv(pnorm(q), "q")
Deriv(dnorm(x, mean), "mean")
Deriv(asin(x), "x")
Deriv(acos(x), "x")
Deriv(atan(x), "x")
Deriv(gamma(x), "x")
Deriv(lgamma(x), "x")
Deriv(digamma(x), "x")
Deriv(trigamma(x), "x")
Deriv(psigamma(x, deriv = 5), "x")
Deriv(x*y, "x")
Deriv(x/y, "x")
Deriv(x^y, "x")
Deriv((x), "x")
Deriv(+x, "x")
Deriv(-x, "x")
Deriv(abs(x), "x")
Deriv(sign(x), "x")

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
