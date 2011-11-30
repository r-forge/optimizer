\name{grnd}
\alias{grnd}
\encoding{UTF-8}
\title{A reorganization of the call to numDeriv grad() function.}
\concept{minimization}
\concept{maximization}
\description{
        Provides a wrapper around user functions for nonlinear optimization
	to try to control for inadmissible arguments to user objective, gradient
	or hessian functions, as well as provide for scaling and maximization.
}
\usage{
       grnd(par, userfn, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{userfn}{A user-supplied function }
 \item{...}{Other data needed to evaluate the user function.}
}
\details{
  None
}
\value{
  \code{grnd} returns an approximation to the gradient of the function fnuser
}
\examples{
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	gg
}

x0<-c(1,2,3,4)
x1<-c(1,1,1,1)
psc<-c(1,2,3,4)
ps1<-c(1,1,1,1)

g0<-genrose.f(x0)
g0
g1<-genrose.f(x1)
g1

tryf01<-grnd(x0, genrose.f, gs=10)
tryf01
cat("Analytic:\n")
g01<-genrose.g(x0,gs=10)
g01
cat("\n\n")
tryf11<-grnd(x1, genrose.f, gs=10)
tryf11
cat("Analytic:\n")
g11<-genrose.g(x1,gs=10)
g11
}
\keyword{nonlinear}
\keyword{optimize}

