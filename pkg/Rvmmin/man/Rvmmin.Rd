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

  This code is
  entirely in R to allow users to explore and understand the method. It also
  allows bounds (or box) constraints and masks (equality constraints) to be
  imposed on parameters. 
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
   function. Note the option of using numerical gradients from package
   \code{numDeriv}.
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
    The source code \code{Rvmmin} for R is still a work in progress, so users should watch
    the console output.

   The \code{control} argument is a list.
   \describe{
   \item{maxit}{A limit on the number of iterations (default 500). Note that this is 
      used to compute a quantity \code{maxfeval}<-round(sqrt(n+1)*maxit) where n is the
      number of parameters to be minimized.}
   \item{maximize}{Set TRUE to maximize the function (default FALSE).}
   \item{trace}{Set 0 (default) for no output, >0 for trace output
      (larger values imply more output).}
   \item{eps}{Tolerance used to calculate numerical gradients. Default is 1.0E-7. See 
      source code for \code{Rvmmin} for details of application.}
   \item{usenumDeriv}{Default is FALSE. Set TRUE to use numDeriv for numerical gradient
      approximations.}
  }
}

\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{fvalue}{The value of the objective at the best set of parameters found.}
  \item{counts}{A vector of two integers giving the number of function and gradient evaluations.}
  \item{convergence}{An integer indicating the situation on termination of the function. \code{0}
   indicates that the method believes it has succeeded. Other values:
   \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
   }
  }
  \item{message}{A description of the situation on termination of the function.}
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
print(ansrosenbrock) # use print to allow copy to separate file that 
#    can be called using source()
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

print(ansbt)

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
cat("genrosea uses analytic gradient\n")
print(genrosea)
cat("genrosenl uses standard numerical gradient\n")
print(genrosenl)
cat("genrosenx uses numDeriv gradient\n")
print(genrosenx)


maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}


negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

cat("test that maximize=TRUE works correctly\n")

n<-6
xx<-rep(1,n)
ansmax<-Rvmmin(xx,maxfn, control=list(maximize=TRUE,trace=1))
print(ansmax)

cat("using the negmax function should give same parameters\n")
ansnegmax<-Rvmmin(xx,negmaxfn, control=list(trace=1))
print(ansnegmax)


#####################
cat("test bounds and masks\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds1<-Rvmmin(startx,genrose.f,lower=lo,upper=up) 
print(grbds1)

cat("test lower bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
grbds2<-Rvmmin(startx,genrose.f,lower=lo) 
print(grbds2)

cat("test lower bound single value only\n")
nn<-4
startx<-rep(pi,nn)
lo<-2
up<-rep(10,nn)
grbds3<-Rvmmin(startx,genrose.f,lower=lo) 
print(grbds3)

cat("test upper bound only\n")
nn<-4
startx<-rep(pi,nn)
lo<-rep(2,nn)
up<-rep(10,nn)
grbds4<-Rvmmin(startx,genrose.f,upper=up) 
print(grbds4)

cat("test upper bound single value only\n")
nn<-4
startx<-rep(pi,nn)
grbds5<-Rvmmin(startx,genrose.f,upper=10) 
print(grbds5)



cat("test masks only\n")
nn<-6
bd<-c(1,1,0,0,1,1)
startx<-rep(pi,nn)
grbds6<-Rvmmin(startx,genrose.f,bdmsk=bd) 
print(grbds6)

cat("test upper bound on first two elements only\n")
nn<-4
startx<-rep(pi,nn)
upper<-c(10,8, Inf, Inf)
grbds7<-Rvmmin(startx,genrose.f,upper=upper) 
print(grbds7)


cat("test lower bound on first two elements only\n")
nn<-4
startx<-rep(0,nn)
lower<-c(0,1.1, -Inf, -Inf)
grbds8<-Rvmmin(startx,genrose.f,genrose.g,lower=lower, control=list(maxit=2000)) 
print(grbds8)

cat("test n=1 problem using simple squares of parameter\n")

sqtst<-function(xx) {
   res<-sum((xx-2)*(xx-2))
}

nn<-1
startx<-rep(0,nn)
onepar<-Rvmmin(startx,sqtst, control=list(trace=1)) 
print(onepar)



}

\keyword{nonlinear}
\keyword{optimize}

