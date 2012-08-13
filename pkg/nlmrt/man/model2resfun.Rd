\name{model2resfun}
\alias{model2resfun}
\title{Generate a residual function from a nonlinear model expression
and a vector of named parameters.
}
\description{
   Given a nonlinear model expressed as an expression of the form
         lhs ~ formula_for_rhs
   and a start vector where parameters used in the model formula are named,
   attempts to build the the \code{R} function for the residuals of the model.
   As a side effect, a text file with the program code is generated.
}
\usage{
   model2resfun(modelformula, pvec, funname="myres", filename=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelformula}{
     This is a modeling formula of the form (as in \code{nls})
           lhsvar ~ rhsexpression
     for example,
           y ~ b1/(1+b2*exp(-b3*tt))
     You may also give this as a string.
}
  \item{pvec}{
     A named parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     WARNING: the parameters in the output function will be used in the 
     order presented in this vector. Names are NOT respected in the 
     output function. 
}
  \item{funname}{
     The (optional) name for the function that is generated in the file
     named in the next argument. The default name is 'myres'.
} 
  \item{filename}{
     The (optional) name of a file that is written containing the (text) program code
     for the function. If NULL, no file is written.
} 
}
\details{
    None.
}
\value{
  An \code{R} function object that computes the gradient of the sum of 
  squared residuals of a nonlinear model at a set of parameters. 
}
\references{
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications
}
\author{
John C Nash <nashjc@uottawa.ca>
}
\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\examples{
  cat("See also examples in nlmrt-package.Rd\n")
    # a test
  y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 
      50.156, 62.948, 75.995, 91.972)  # for testing
  tt <- seq_along(y)  # for testing
  # NOTE: use of t gives confusion with t() in R CMD check,
  # but not in direct use with source() 120429
  f <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt))
  p <- c(b1 = 1, b2 = 1, b3 = 1)
  myres <- model2resfun(f, p)
  cat("myres:\n")
  print(myres)
  ans <- myres(p, tt = tt, y = y)  
  cat("ans:")  
  print(ans)
  cat("ss ( =? 23520.58):", as.numeric(crossprod(ans)),"\n")
  bnew <- c(b1 = 200, b2 = 50, b3 = 0.3)
  # anew<-myres(prm=bnew, t=t, y=y) # Note issue with t vs
  # t()
  anew <- eval(myres(prm = bnew, tt = tt, y = y))
  cat("anew:")
  print(anew)
  cat("ss ( =? 158.2324):", as.numeric(crossprod(anew)),"\n")
  cat("Test with vector of un-named parameters\n")
  bthree <- c(100, 40, 0.1)
  athree <- try(myres(prm = bthree, tt = tt, y = y))
  print(athree)
  cat("ss ( =? 19536.65):", as.numeric(crossprod(athree)),"\n")
}
\keyword{ nonlinear least squares }

