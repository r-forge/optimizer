\name{Simplify}
\alias{newSimplification}
\alias{Simplify}
\alias{sysSimplifications}
\alias{isFALSE}
\alias{isZERO}
\alias{isONE}
\alias{isMINUSONE}
\title{
Functions to simplify expressions.
}
\description{
\code{Simplify} simplifies expressions according to rules specified
by \code{newSimplification}.
}
\usage{
Simplify(expr, simpEnv = sysSimplifications, verbose = FALSE)
newSimplification(expr, test, simplification, do_eval = FALSE, simpEnv = sysSimplifications)
sysSimplifications
isFALSE(x)
isZERO(x)
isONE(x)
isMINUSONE(x)
}
\arguments{
  \item{expr}{
An expression to simplify; this should be a language object.
}
  \item{simpEnv}{
An environment holding the registered simplifications.
}  
  \item{verbose}{
If \code{TRUE}, diagnostics are printed as simplifications are recognized.
}
  \item{test}{
An expression giving a test to apply to decide whether this simplification
applies.
}
  \item{simplification}{
The new expression to use to replace the original one.
}
  \item{do_eval}{
Whether to evaluate the new expression (to give an even newer
expression) to use as the simplification.
}

  \item{x}{
An expression to test.
}
}
\details{
The \code{Simplify} function uses simple rules to simplify
expressions.  The simplification is aimed at the needs of this
package, so the built-in simplificatinos assume that variables and expressions 
have finite real values. For example, \code{0*expr} will simplify to \code{0}
regardless of the value of \code{expr}.

\code{newSimplification} adds a new simplification pattern to the
registered collection.  The tests are applied to function calls with
the same function and number of parameters, in order as specified.
Users may specify their own environment (perhaps parented by
\code{sysSimplifications}) to hold rules if they wish to override the
standard rules.

The \code{isFALSE}, \code{isZERO}, \code{isONE}, and \code{isMINUSONE}
functions are simple functions to test whether expressions are simple constants,
similar to \code{\link{isTRUE}}.

To handle functions which act differently depending on which arguments
are present, \code{Simplify} will simplify the expression \code{missing(.MissingVal)}
to \code{TRUE}.  This is used in the definition of the derivative
for \code{x - y}, where the unary minus is seen as a missing \code{y}
value.
}
\value{
\code{Simplify} returns a simplification of the expression, assuming
that variables and functions take real values.

The \code{newSimplification} function is called for the side
effect of recording a new simplification rule.

If \code{expr} or \code{deriv} is missing, \code{newSimplification()} will
report on the currently saved simplifications in \code{simpEnv}.
}
\author{
Duncan Murdoch
}
\seealso{
\code{\link{Deriv}}, which makes use of \code{Simplify}.
}
\examples{
Simplify(quote(a + 0))
Simplify(quote(exp(1)), verbose = TRUE)

Simplify(quote(sqrt(a + b)))  # standard rule
myrules <- new.env(parent = sysSimplifications)
newSimplification(sqrt(a), TRUE, a^0.5, simpEnv = myrules)
Simplify(quote(sqrt(a + b)), simpEnv = myrules)
}
\keyword{ math }
