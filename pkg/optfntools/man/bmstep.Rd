\name{bmstep}
\alias{bmstep}
\encoding{UTF-8}
\title{Compute the maximum step along a search direction.}
\description{Nonlinear optimization problems often have explicit or implicit upper and
lower bounds on the parameters of the function to be miminized or maximized. These are 
called bounds or box constraints. Some of the parameters may be fixed for a given problem
or for a temporary trial. These fixed, or masked, paramters are held at one value during 
a specific 'run' of the optimization.

The bmstep() function computes the maximum step possible (which could be infinite) 
along a particular search direction from current parameters to bounds.

}
\usage{
   bmstep(par, srchdirn, lower=NULL, upper=NULL, bdmsk=NULL, trace=0)
}
\arguments{
 \item{par}{A numeric vector of starting values of the optimization function parameters.}
 \item{srchdirn}{A numeric vector giving the search direction.}
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
     The bmstep function will compute and return (as a double or Inf) the 
     maximum step to the bounds.
}

\value{
  A double precision value or Inf giving the 
     maximum step to the bounds.
}
\examples{
#####################
# tbmstep.R
#####################
# rm(list=ls())
# source("/home/john/R-optimtest/work/optfntools/R/bmstep.R")
cat("Test bmstep() function \n")
# test bmstep ideas
tdir<-c(1,0,-1) # search direction
up<-c(Inf, 10, 10) # upper bound
lo<-c(-10, -10, -Inf) # lower bound
bvec<-c(1,1,1) # current parameters
bm<-c(1,1,1)

dall<-data.frame(bvec, tdir,up, lo)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=NULL, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

cat("\n\n change direction of search\n")
tdir<- -c(1,0,-1) # search direction
dall<-data.frame(bvec, tdir,up, lo, bm)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=bm, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

cat("\n Now mask first parameter\n")
bm<-c(0,1,1)
dall<-data.frame(bvec, tdir,up, lo, bm)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=bm, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

cat("\n Now mask first parameter and use to original search direction\n")
bm<-c(0,1,1)
tdir<-c(1,0,-1) # search direction
dall<-data.frame(bvec, tdir,up, lo, bm)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=bm, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

cat("\nMask last parameter and use original search direction\n")
bm<-c(1,1,0)
tdir<-c(1,0,-1) # search direction
dall<-data.frame(bvec, tdir,up, lo, bm)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=bm, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

cat("\n Now mask last parameter and neg original search direction\n")
bm<-c(1,1,0)
tdir<- -c(1,0,-1) # search direction
dall<-data.frame(bvec, tdir,up, lo, bm)

ms<-bmstep(bvec, srchdirn=tdir, lower=lo, upper=up, bdmsk=bm, trace=2)
print(dall)
cat("maxst=",ms,"\n\n")

}
\keyword{nonlinear}
\keyword{optimize}
\keyword{upper}
\keyword{lower}
\keyword{bound}
\keyword{mask}
