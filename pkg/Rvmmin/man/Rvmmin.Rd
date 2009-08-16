\name{Rvmmin}
\alias{Rvmmin}
\encoding{UTF-8}
\title{Variable metric nonlinear function minimization with bounds constraints}
\description{An R implementation of a nonlinear variable metric method for minimization
of nonlinear functions possibly subject to bounds (box) constraints and masks 
(fixed parameters). The algorithm is based on Nash (1979) Algorithm 21 for main structure,
which is itself drawn from Fletcher's (1970) variable metric code. This is also the basis
of optim() method 'BFGS' which, however, does not deal with bounds or masks. In the 
present method, an approximation to the inverse Hessian (B) is used to generate a 
search direction t = - B &*& g, a simple backtracking line search is used until an
acceptable point is found, and the matrix B is updated using a BFGS formula. If no
acceptable point can be found, we reset B to the identity i.e., the search direction
becomes the negative gradient. If the search along the negative gradient is unsuccessful,
the method terminates. 
}
\usage{
   Rvmmin(par, fn, gr, lower, upper, bdmsk, control = list(), \dots)
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
   is discouraged, as the method seems to work better with accurate
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
    Successful completion.
    The codes \code{Rvmmin} for R are still a work in progress, so users should watch
    the console output.
  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
  There are other items to be added.??
}
\references{
  Fletcher, R (1970) ...

  Nash, J C (1979, 1990) ...

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
ansrosenbrock <- Rvmmin(fn=fr,par=c(1,2))
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
ansbt<-Rvmmin(xx, bt.f, bt.g, lower, upper, bdmsk, control=list(trace=1))

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
genrosea<-Rvmmin(xx,genrose.f, genrose.g)
genrosenl<-Rvmmin(xx,genrose.f) # use local numerical gradient
genrosenx<-Rvmmin(xx,genrose.f,control=list(usenumDeriv=TRUE)) # use numDeriv for gradients
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
ansmax<-Rvmmin(xx,maxtest.f, maxtest.g,control=list(maximize=TRUE,trace=1))
ansmaxn<-Rvmmin(xx,maxtest.f, control=list(maximize=TRUE,trace=1))


}
\keyword{nonlinear}
\keyword{optimize}

