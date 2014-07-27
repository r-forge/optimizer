\name{Deriv}
\alias{Deriv}
\alias{newDeriv}
\alias{sysDerivs}
\title{
Functions to take symbolic derivatives.
}
\description{
\code{Deriv} takes a symbolic derivative of the first argument.  
}
\usage{
Deriv(expr, name, envir = sysDerivs, do_substitute = TRUE)
newDeriv(expr, deriv, envir = sysDerivs) 
sysDerivs
}
\arguments{
  \item{expr}{
An expression vector or language object.
}
  \item{name}{
The name of the variable with respect to which the derivative will be computed.
}
  \item{envir}{
The environment in which derivatives are stored.
}
  \item{do_substitute}{
If \code{TRUE}, use \code{\link{substitute}} to get the expression passed as
\code{expr}, otherwise evaluate it.
}
  \item{deriv}{
An expression giving the derivative of the function call in \code{expr}.
}
}
\details{
The \code{Deriv} function computes symbolic derivatives of expressions
or language objects. The recognized derivatives are stored in
\code{envir}; the default \code{sysDerivs} contains expressions for
all of the derivatives recognized by \code{\link{deriv}}, but in
addition allows differentiation with respect to any parameter
where it makes sense.  It also allows the derivative of \code{abs}
and \code{sign}, using an arbitrary choice of 0 at the discontinuities.

The \code{newDeriv} function can be used to define a new derivative.
The \code{expr} argument should match the header of the function as a
call to it (e.g. as in the help pages), and the \code{deriv} argument
should be an expression giving the derivative, including calls to
\code{D(arg)}, which will not be evaluated, but will be substituted
with partial derivatives of that argument with respect to \code{name}.
See the examples below.

If \code{expr} or \code{deriv} is missing in a call to
\code{newDeriv()}, it will return the currently saved derivative
record from \code{envir}.  If \code{name} is missing in a call to
\code{Deriv} with a function call, it will print a message describing
the derivative formula and return \code{NULL}.

To handle functions which act differently if a parameter is
missing, code the default value of that parameter to \code{.MissingVal},
and give a derivative that is conditional on \code{missing()}
applied to that parameter.  See the derivative of \code{"-"}
for an example.
}
\value{
If \code{expr} is an expression vector, \code{Deriv} and \code{Simplify}
return expression vectors containing the response.

\code{Deriv} returns the symbolic derivative of the expression.

\code{newDeriv} with \code{expr} and \code{deriv} specified is
called for the side effect of recording the derivative in \code{envir}.
If \code{expr} is missing, it will return the list of names of functions
for which derivatives are recorded.  If \code{deriv} is missing, it
will return its record for the specified function.
}
\author{
Duncan Murdoch
}
\seealso{
\code{\link{deriv}}, \code{\link{Simplify}}
}
\examples{
newDeriv()
newDeriv(sin(x))
Deriv(sin(x+y), "x")

f <- function(x) x^2
newDeriv(f(x), 2*x*D(x))
Deriv(f(abs(x)), "x")

Deriv(pnorm(x, sd=2, log = TRUE), "x")
}
\keyword{math}
\keyword{nonlinear}
