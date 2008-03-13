\name{spg}
\alias{spg}
\title{Spectral Projected Gradient}
\description{Spectral projected gradient method for
  	  box-constrained minimization.}
\usage{
     spg(par, fn, grad=NULL, lower=-Inf, upper=Inf,  
	grad.method="simple", control=list(), ...)
 }
\arguments{
    \item{fn}{a function with a scalar real result (see details).}
    \item{par}{DOCUMENTa real scalar or vector argument to func, indicating the 
      point(s) at which the dfsaneient is to be calculated.}
    \item{grad}{DOCUMENT} 
    \item{lower}{DOCUMENT} 
    \item{upper}{DOCUMENT} 
    \item{grad.method}{DOCUMENT} 
    \item{control}{A list of control parameters. See *Details*.}
    \item{...}{an additional arguments passed to \code{func}.}
}
\value{FOR EXAMPLEA real scalar or vector of the approximated spgient(s).}

\details{
 R translation (with minor modifications):  Ravi Varadhan, Johns Hopkins University, February 23, 2008.
 Incorporates box constraints

Argument \code{control} is list specifing any changes to default values of
algorithm control parameters. The default values are
    \code{M=10},
    \code{maxit=2500} 
    \code{tol=1e-08} 
    \code{gtol=1e-04} 
    \code{maxfeval=10000} 
    \code{trace=FALSE} DOCUMENT WHAT THESE DO

FOR EXAMPLE
   The function \code{spg} calculates a numerical approximation of the 
   first derivative of \code{func} at the point \code{x}. Any additional 
   arguments in \dots are also passed to \code{func}, but the spgient is not
   calculated with respect to these additional arguments.
   It is assumed \code{func}  is a scalar value function. If a vector \code{x} 
   produces a scalar
   result then \code{spg} returns the numerical approximation of the spgient
   at the point \code{x} (which has the same length as \code{x}).
   If a vector \code{x} produces a vector result then the result must have the
   same length as \code{x}, and it is assumed that this corresponds to applying
   the function to each of its arguments (for example, \code{sin(x)}). 
   In this case \code{spg} returns the
   spgient at each of the points in \code{x} (which also has the same length 
   as \code{x} - so be careful). An alternative for vector valued functions is
   provided by.
  
   If method is "simple", the calculation is done using a simple epsilon
   difference. For this case, only the element \code{eps} of \code{methods.args}
   is used.

   If method is "Richardson", the calculation
   is done by Richardson's extrapolation (see e.g. Linfield and Penny, 1989,
   or Fornberg and Sloan, 1994.)
   This method should be used if accuracy, as opposed to speed, is important.
   For this case, 
   \code{methods.args=list(eps=1e-4, d=0.01, r=4, show.details=FALSE)} 
   are used.
    \code{d} gives the fraction of \code{x} to use for the initial numerical 
      approximation. The default means the initial approximation uses
      \code{0.0001 * x}.
    \code{eps} is used instead of \code{d} for elements of \code{x} which are 
       zero.
    \code{r} gives the number of Richardson improvement iterations (repetitions
       with successly smaller \code{d}. The default \code{4} general provides 
       good results, but this can be increased to \code{6} for improved
       accuracy at the cost of more evaluations.
    \code{v} gives the reduction factor.
    \code{show.details} is a logical indicating if detailed calculations should 
        be shown.
   
  The general approach in the Richardson method is to iterate for \code{r} 
  iterations from initial 
  values for interval value \code{d},  using reduced factor \code{v}.
  The the first order approximation to the derivative with respect 
  to \eqn{x_{i}}{x_{i}} is

      \deqn{f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}{%
	    f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}
       
  This is repeated \code{r} times  with successively smaller \code{d}  and 
	  then Richardson extraplolation. is applied.

}
\references{ 
    Non-monotone spectral projected-gradient method for minimization
    Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected 
    gradient methods on convex sets, \emph{SIAM J Optimization}, 10, 1196-1211.
    
    Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for 
    convex-constrained optimization, 
    \emph{ACM Transactions on Mathematical Software}.
}
\seealso{
  \code{\link{sane}},
  \code{\link{dfsane}},
  \code{\link{nlsolve}},
  \code{\link[stats]{numericDeriv}}
   }
\examples{
trig.f <- function(x){
   f <- length(x) - sum(cos(x)) + (1:length(x))*(1 - cos(x)) - sin(x) 
   sum(f*f)
   }

ans.spg <- spg(par=runif(100), fn=trig.f, maxit=2500)
}
\keyword{multivariate}
