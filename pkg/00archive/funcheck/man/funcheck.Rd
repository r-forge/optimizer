\name{funcheck}


\alias{funcheck}
\alias{funtest}


\title{Run tests, where possible, on functions to be optimized}


\description{

   There are two functions provided, \code{funcheck} and \code{funtest}. 

   \code{funcheck} works on a source code file containing several functions with
   a common root to their name. Thus, a file 'hobbs.R' may contain code for several
   functions:
   \describe{
   \item{hobbs.f}{the objective function to be minimized or maximized. This takes as 
     arguments a vector of parameters and local variables specified via the "..." 
     inputs.}
   \item{hobbs.g}{the gradient of that function}
   \item{hobbs.h}{the hessian of that funciton}
   \item{hobbs.res}{if hobbs.f is formed as a sum of squares, then this function returns
       a vector of residuals such that the sum of their squares is the function value}
   \item{hobbs.jac}{to compute the Jacobian of the residuals above}
   \item{hobbs.rsd}{to compute the second derivatives of the residuals above}
   \item{hobbs.setup}{to provide starting values and bounds for the functions above}
   }

   For \code{funtest}, the functions are specified explicitly.

}


\usage{
funcheck(xpar, fname, lower=NULL, upper=NULL, cctrl=list(trace=1), ... )
funtest(xpar, fn, gr=NULL, hess=NULL, res=NULL, jac=NULL, rsd=NULL, lower=NULL, upper=NULL, cctrl=list(trace=1), ... )
}


\arguments{
     
    \item{xpar}{
        a matrix with rows each being a vector giving the numerical values of an initial 
        point where the function should be tested.
        }
    \item{fname}{
	the base name (as a character string) of a file of type .R that contains functions fname.f 
        to compute the scalar 
        objective function value at a given set of parameters in the vector x, where x is of length n, 
        with (optionally) fname.g for the gradient, fname.res for residuals in the case that the 
	scalar objective function is a sum of m squared elements, fname.jac for the Jacobian of these
	residuals, fname.setup to set up the xset, as well as lower and upper bounds and the bounds
	and masks indicator (integer) vector bdmsk, and fname.h for the Hessian. 
        }
    \item{fn}{explicit name for the objective function}
    \item{gr}{explicit name for the gradient function}
    \item{hess}{explicit name for the Hessian function}
    \item{res}{explicit name for the residual function}
    \item{jac}{explicit name for the Jacobian function}
    \item{rsd}{explicit name for the residual second derivative function}
    \item{lower}{
	a matrix whose rows are sets of lower bounds on the parameters to the function
        }
    \item{upper}{
	a matrix whose rows are sets of upper bounds on the parameters to the function
        }
    \item{cctrl}{
        list of control parameters.\cr
	There are two elements currently: \cr
		trace: set >0 to provide output from funcheck to the console, 0 otherwise
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }

}


\details{
\tabular{ll}{
Package: \tab funcheck\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
    \bold{Algorithm}

	\bold{funcheck} assumes that the user supplies functions in the file fname.R
	that is accessible from the global environment.

        \bold{funtest} assumes that the functions called are already in the current
        environment.

}


\value{A list of lists, one list for each set of parameters (each row) of xset. Each sublist has
the following components:

\item{fname}{The name of the test function file.}

\item{ipoint}{An index counter for the test points.}

\item{par}{The numerical values of the n elements of the supplied parameters..} 

\item{lower}{Lower bounds for the parameters in this case. }

\item{upper}{Lower bounds for the parameters in this case. }

\item{bstate}{A character vector telling whether each parameter is Out of Bounds Low, In Bounds, or Out
    of Bounds High. }

\item{grd}{The max. abs. relative difference with denominator augmented by 10. for the gradient.}

\item{Jrd}{The max. abs. relative difference with denominator augmented by 10. for the Jacobian.}

\item{jsdrd}{The max. abs. relative difference with denominator augmented by 10. for the Jacobian 
     second derivatives.}

\item{Hrd}{The max. abs. relative difference with denominator augmented by 10. for the Hessian.}

}


\note{
	These versions of funcheck and funtest are not optimized.

The following code must be stored in file tryfun.R and available to the tests. ??how?? ??where??


}



\author{

    John C. Nash
}


\examples{
## the file tryfun.R must be available in the current directory
#don't execute --   ans<-funcheck(xpar = rep(10.0,10), fname="tryfun") 
# The above example is in the tests

}

\keyword{optimize}

