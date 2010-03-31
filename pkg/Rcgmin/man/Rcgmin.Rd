\name{Rcgmin}
\alias{Rcgmin}
\encoding{UTF-8}
\title{An R implementation of a nonlinear conjugate gradient algorithm with the Dai / Yuan
   update and restart. Based on Nash (1979) Algorithm 22 for main structure.}
\description{
  The purpose of \code{Rcgmin} is to minimize an unconstrained or bounds
  (box) and mask constrained function 
  of many parameters by a nonlinear conjugate gradients method. 
}
\usage{
   Rcgmin(par, fn, gr, lower, upper, bdmsk, control = list(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{gr}{A function that returns the gradient of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. This function 
   returns the gradient as a numeric vector.

   If gr is NOT supplied, a numerical approximation will be used. This 
   is STRONGLY discouraged, as the CG methods seem to require accurate
   gradients. Furthermore, the computation of numerical gradients is
   many times slower -- of the order of 1000s -- than using an analytic
   function.
}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization.}
 \item{control}{
    An optional list of control settings.  
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
   Functions \code{fn} must return a numeric value.

   The \code{control} argument is a list.
   \describe{
   \item{maxit}{A limit on the number of iterations (default 500). Note that this is 
      used to compute a quantity \code{maxfeval}<-round(sqrt(n+1)*maxit) where n is the
      number of parameters to be minimized.}
   \item{maximize}{Set TRUE to maximize the function (default FALSE).}
   \item{trace}{Set 0 (default) for no output, >0 for trace output
      (larger values imply more output).}
    \item{}{The source code \code{Rcgmin} for R is likely to remain a work in progress for some time, 
    so users should watch the console output.}
  }
}
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of the objective at the best set of parameters found.}
  \item{counts}{A two-element integer vector giving the number of calls to
          'fn' and 'gr' respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to 'fn'
          to compute a finite-difference approximation to the gradient.}
  \item{convergence}{An integer code. 
	 '0' indicates successful convergence.
         '1' indicates that the function evaluation count 'maxfeval' was reached.
         '2' indicates initial point is infeasible.}
  \item{message}{A character string giving any additional information returned
          by the optimizer, or 'NULL'.}
}
\references{
To be added.
}
\seealso{\code{\link{optim}}}
\examples{
#####################
   ## Rosenbrock Banana function
fr <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
ansrosenbrock <- Rcgmin(fn=fr,par=c(1,2))
ansrosenbrock
#####################
# Simple bounds and masks test
bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

n<-10
xx<-rep(0,n)
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
   lower[i]<-1.0*(i-1)*(n-1)/n
   upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)
ansbt<-Rcgmin(xx, bt.f, bt.g, lower, upper, bdmsk, control=list(trace=1))

#####################
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}
genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
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

# analytic gradient test
xx<-rep(pi,10)
lower<-NULL
upper<-NULL
bdmsk<-NULL
genrosea<-Rcgmin(xx,genrose.f, genrose.g)
genrosenl<-Rcgmin(xx,genrose.f) # use local numerical gradient
genrosenx<-Rcgmin(xx,genrose.f,control=list(usenumDeriv=TRUE)) # use numDeriv for gradients
genrosea
genrosenl
genrosenx

# How about some maximization tests??
# maxtest function -- a function to maximize
# ref??


maxtest.f <- function(x) {
 n <- length(x)
 res<-(seq(3,(n+2))+0.25*rep(1,n)-x)
 ff<-exp(-sum(res*res)/25)
 ff
 }

maxtest.g <- function(x) {
 n <- length(x)
 res<-seq(3,(n+2))+0.25*rep(1,n)-x
 g<- 2*res*exp(-sum(res*res)/25)/25
 g
}

n<-6
xx<-rep(1,n)
ansmax<-Rcgmin(xx,maxtest.f, maxtest.g,control=list(maximize=TRUE,trace=1))
ansmaxn<-Rcgmin(xx,maxtest.f, control=list(maximize=TRUE,trace=1))


}



\keyword{nonlinear}
\keyword{optimize}

