\name{model2jacfun}
\alias{model2jacfun}
\title{Generate a Jacobian matrix function from a nonlinear model expression
and a vector of named parameters.
}
\description{
   Given a nonlinear model expressed as an expression of the form
         lhs ~ formula_for_rhs
   and a start vector where parameters used in the model formula are named,
   attempts to build the the \code{R} function for the Jacobian of the residuals.
   As a side effect, a text file with the program code is generated.
}
\usage{
   model2jacfun(modelformula, pvec, funname="myjac", filename=NULL)
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
     named in the next argument. The default name is 'myjac'.
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
  An \code{R} function object that computes the Jacobian of the nonlinear
  model at a set of parameters.
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
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
    tt<-1:length(y) # for testing
    f<- y ~ b1/(1+b2*exp(-1*b3*tt))
    p<-c(b1=1, b2=1, b3=1)
    myfn<-model2jacfun(f, p)
    cat("myfn:\n")
    print(myfn)
    
    ans<-myfn(p, tt=tt,y=y)
    print(ans)
    bnew<-c(b1=200, b2=50, b3=0.3)
    anew<-myfn(prm=bnew, tt=tt, y=y)
    print(anew)
   
    cat("Test with un-named vector\n") # At 20120424 should fail
    bthree<-c(100,40,.1)
    athree<-try(myfn(prm=bthree, tt=tt, y=y))
    print(athree)
}
\keyword{ nonlinear least squares }

