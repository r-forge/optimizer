% File src/library/stats/man/optimx.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{opmxtra}
\alias{opmxtra}
% \alias{[.optimx}
% \alias{as.data.frame.optimx}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  A function to extend the results dataframe of the \code{opm}
  multi-method optimization function; in particular, to compute and
  report the Kuhn-Karush-Tucker tests.
}
\usage{
opmxtra <- function(ansout, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            control=list(), ...) 
}
\arguments{
 \item{ansout}{a dataframe of results from the \code{opm} function.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.
}
\value{

   If there are \code{npar} parameters, then the result is a dataframe having one row
   for each method for which results are reported, using the method as the row name,
   with columns

   \code{par_1, .., par_npar, value, fevals, gevals, niter, convcode, kkt1, kkt2, xtimes}

  where
  \describe{
  \item{par_1}{ .. }
  \item{par_npar}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{fevals}{The number of calls to \code{fn}.}
  \item{gevals}{The number of calls to \code{gr}. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.}
  \item{niter}{For those methods where it is reported, the number of ``iterations''. See
    the documentation or code for particular methods for the meaning of such counts.}
  \item{convcode}{An integer code. \code{0} indicates successful
    convergence. Various methods may or may not return sufficient information
	to allow all the codes to be specified. An incomplete list of codes includes
    \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
      \item{\code{10}}{indicates degeneracy of the Nelder--Mead simplex.}
      \item{\code{51}}{indicates a warning from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{52}}{indicates an error from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
    }
  }
  \item{kkt1}{A logical value returned TRUE if the solution reported has a ``small'' gradient.}
  \item{kkt2}{A logical value returned TRUE if the solution reported appears to have a 
  positive-definite Hessian.}
  \item{xtimes}{The reported execution time of the calculations for the particular method.}
  }


The attribute "details" to the returned answer object contains information,
if computed, on the gradient (\code{ngatend}) and Hessian matrix (\code{nhatend}) 
at the supposed optimum, along with the eigenvalues of the Hessian (\code{hev}), 
as well as the \code{message}, if any, returned by the computation for each \code{method},
which is included for each row of the \code{details}. 
If the returned object from optimx() is \code{ans}, this is accessed 
via the construct
    \code{attr(ans, "details")}

This object is a  matrix based on a list so that if ans is the output of optimx
then attr(ans, "details")[1, ] gives the first row and 
attr(ans,"details")["Nelder-Mead", ] gives the Nelder-Mead row. There is 
one row for each method that has been successful 
or that has been forcibly saved by save.failures=TRUE. 

There are also attributes
   \describe{
   \item{maximize}{to indicate we have been maximizing the objective}
   \item{npar}{to provide the number of parameters, thereby facilitating easy
        extraction of the parameters from the results data frame}
   \item{follow.on}{to indicate that the results have been computed sequentially,
        using the order provided by the user, with the best parameters from one
        method used to start the next. There is an example (\code{ans9}) in 
        the script \code{ox.R} in the demo directory of the package.}
   }
}
\source{

See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

}
\references{

 See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

 Nash JC (2014). On Best Practice Optimization Methods in R., 
        \emph{Journal of Statistical Software}, 60(2), 1-14.,
        URL http://www.jstatsoft.org/v60/i02/.

}

\examples{
require(graphics)
cat("Note possible demo(ox) for extended examples\n")


## Show multiple outputs of optimx using all.methods
# genrose function code
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
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

startx<-4*seq(1:10)/3.
ans8<-opm(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
   control=list(trace=0), gs=10)
ans8

ans8x <- opmxtra(ans8,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
   control=list(save.failures=TRUE, trace=0), gs=10)


ans8[, "gevals"]
ans8["spg", ]
summary(ans8, par.select = 1:3)
summary(ans8, order = value)[1, ] # show best value
head(summary(ans8, order = value)) # best few
## head(summary(ans8, order = "value")) # best few -- alternative syntax

## order by value.  Within those values the same to 3 decimals order by fevals.
## summary(ans8, order = list(round(value, 3), fevals), par.select = FALSE)
summary(ans8, order = "list(round(value, 3), fevals)", par.select = FALSE)

## summary(ans8, order = rownames, par.select = FALSE) # order by method name
summary(ans8, order = "rownames", par.select = FALSE) # same

summary(ans8, order = NULL, par.select = FALSE) # use input order
## summary(ans8, par.select = FALSE) # same

}
\keyword{nonlinear}
\keyword{optimize}
