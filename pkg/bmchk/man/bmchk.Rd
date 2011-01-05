\name{bmchk}
\alias{bmchk}
\encoding{UTF-8}
\title{Check bounds and masks for parameter constraints used in nonlinear optimization}
\description{Nonlinear optimization problems often have explicit or implicit upper and
lower bounds on the parameters of the function to be miminized or maximized. These are 
called bounds or box constraints. Some of the parameters may be fixed for a given problem
or for a temporary trial. These fixed, or masked, paramters are held at one value during 
a specific 'run' of the optimization.

It is possible that the bounds are inadmissible, that is, that at least one lower bound
exceeds an upper bound. In this case we set the flag \code{admissible} to FALSE.

Parameters that are outside the bounds are moved to the nearest bound and the flag
\code{parchanged} is set TRUE. However, we DO NOT change masked parameters, and they
may be outside the bounds. This is an implementation choice, since it may be useful
to test objective functions at point outside the bounds.

The package bmchk is essentially a test of the R function bmchk(), which is likely to be 
incorporated within optimization codes.
}
\usage{
   bmchk(par, lower=NULL, upper=NULL, bdmsk=NULL, trace=0)   
}
\arguments{
 \item{par}{A numeric vector of starting values of the optimization function parameters.}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization. Partly for historical reasons, we use the 
     same array during the progress of optimization as an indicator that a 
     parameter is at a lower bound (bdmsk element set to -3) or upper bound (-1).}
 \item{trace}{An integer that controls whether diagnostic information is displayed.
              A positive value displays information, 0 (default) does not.}
}
\details{
     The bmchk function will check that the bounds exist and are admissible, 
     that is, that there are no lower bounds that exceed upper bounds. 

     NOTE: Free paramters outside bounds are adjusted to the nearest bound.
     We then set parchanged = TRUE which implies the original parameters
     were infeasible.

     There is a check if lower and upper bounds are very close together, in 
     which case a mask is imposed and maskadded is set TRUE. NOTE: it is 
     generally a VERY BAD IDEA to have bounds close together in optimization,
     but here we use a tolerance based on the double precision machine 
     epsilon. Thus it is not a good idea to rely on bmchk() to test if 
     bounds constraints are well-posed.
 
}

\value{
  A list with components:
  \item{bvec}{The vector of parameters, possibly adjusted for bounds. Parameters 
         outside bounds are adjusted to the nearest bound.}
  \item{nolower}{TRUE if no lower bounds, FALSE otherwise}
  \item{noupper}{TRUE if no upper bounds, FALSE otherwise}
  \item{bounds}{TRUE if there are any bounds, FALSE otherwise}
  \item{admissible}{TRUE if bounds are admissible, FALSE otherwise
         This means no lower bound exceeds an upper bound. That is the bounds 
         themselves are sensible. This condition has nothing to do with the 
         starting parameters.}
  \item{maskadded}{TRUE when a mask has been added because bounds are very close
        or equal, FALSE otherwise. See the code for the implementation.}
  \item{parchanged}{TRUE if parameters are changed by bounds, FALSE otherswise.
        Note that parchanged = TRUE implies the input parameter values were infeasible, 
        that is, violated the bounds constraints.}
}
\examples{
#####################
# rm(list=ls())
cat("Test bmchk() function \n")
require(bmchk) # load the sub-program

# try free parameters and bounds
prm<-rep(pi,4)
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)

a1<-bmchk(prm, lo, up, bm, trace=1)
cat("a1\n")
print(a1)
cat("\n\n")

# set one parameter low
prm<-rep(pi,4)
prm[3]<-0
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)

a2<-bmchk(prm, lo, up, bm, trace=1)
cat("a2\n")
print(a2)
cat("\n\n")

# set one parameter hi
prm<-rep(pi,4)
prm[2]<-10
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)

a3<-bmchk(prm, lo, up, bm, trace=1)
cat("a3\n")
print(a3)
cat("\n\n")

# try hi and low

prm<-rep(pi,4)
prm[3]<-0
prm[2]<-10
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)

a4<-bmchk(prm, lo, up, bm, trace=1)
cat("a4\n")
print(a4)
cat("\n\n")

# try masked parameters and bounds
prm<-rep(pi,4)
prm[3]<-0
prm[2]<-10
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)
bm[2]<-0
bm[3]<-0

a5<-bmchk(prm, lo, up, bm, trace=1)
cat("a5\n")
print(a5)
cat("\n\n")

# try bounds close together
prm<-rep(pi,4)
lo<-rep(1,4)
up<-rep(5,4)
bm<-rep(1,4)
lo[2]=up[2]-0.001*.Machine$double.eps

a6<-bmchk(prm, lo, up, bm, trace=1)
cat("a6\n")
print(a6)
cat("\n\n")

# try bounds inadmissible
prm<-rep(pi,4)
lo<-rep(2,4)
up<-rep(3,4)
bm<-rep(1,4)
lo[2]<-6


a7<-bmchk(prm, lo, up, bm, trace=1)
cat("a7\n")
print(a7)
cat("\n\n")



}
\keyword{nonlinear}
\keyword{optimize}
\keyword{upper}
\keyword{lower}
\keyword{bound}
\keyword{mask}
