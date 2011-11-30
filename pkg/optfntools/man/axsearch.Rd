\name{axsearch}
\alias{axsearch}
\encoding{UTF-8}
\title{Perform axial search around a supposed minimum and provide diagnostics}
\description{Nonlinear optimization problems often terminate at points in the 
parameter space that are not satisfactory optima. This routine conducts an axial
search, stepping forward and backward along each parameter and computing the objective
function. This allows us to compute the \code{tilt} and \code{radius of curvature} or
\code{roc} along that parameter axis. 

(we will
discuss the problem always in terms of minimization, but the control \code{maximize}
permits maximization. ?? need to expand.)

It is possible that a function value lower than that at the supposed minimum is found.
In this case \code{axsearch} exits immediately with the new function value and
parameters.

}
\usage{
       axsearch(par, fn=NULL, fmin=NULL, lower=NULL, upper=NULL, bdmsk=NULL, trace=0, ...)
}
\arguments{
 \item{par}{A numeric vector of values of the optimization function parameters that are
     at a supposed minimum.}
 \item{fn}{The user objective function}
 \item{fmin}{The value of the objective function at the parameters \code{par}. ?? what if fmin==NULL?}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization. Partly for historical reasons, we use the 
     same array during the progress of optimization as an indicator that a 
     parameter is at a lower bound (bdmsk element set to -3) or upper bound (-1).}
 \item{trace}{If trace>0, then local output is enabled.}
 \item{...}{Extra arguments for the user function.}
}


\details{
None.
}

\value{
  A list with components:
  \item{bestfn}{The lowest (best) function value found (??maximize??) during the axial search, 
        else the original fmin value. (This is actively set in that case.)}
  \item{par}{The vector of parameters at the best function value. }
  \item{details}{A data frame reporting the original parameters, the forward step and backward
       step function values, the size of the step taken for a particular parameter, the tilt and
       the roc (radius of curvature). Some elements will be NA if we find a lower function
       value during the axial search.}
}
\examples{
#####################
sq<-function(x){
   nn<-length(x)
   yy<-1:nn
   f<-sum((yy-x)^2)
}

cat("axial search test at the minimum of a simple function\n")

xx<-1:4
xx
res1<-axsearch(xx, fn=sq)
res1

cat("axial search test NOT at the minimum of a simple function\n")
xy<-rep(pi,3)
res2<-axsearch(xy, fn=sq)
res2

# ?? put in maxfn as a test -- Need to make mods to the program for this.



}
\keyword{nonlinear}
\keyword{optimize}
\keyword{axial}
\keyword{search}
