\name{nlfb}
\alias{nlfb}
\title{Nash variant of Marquardt nonlinear least squares solution via
     qr linear solver.
}
\description{
   Given a nonlinear model expressed as an expression of the form
         lhs ~ formula_for_rhs
   and a start vector where parameters used in the model formula are named,
   attempts to find the minimum of the residual sum of squares using the
   Nash variant (Nash, 1979) of the Marquardt algorithm, where the linear 
   sub-problem is solved by a qr method.
}
\usage{
   nlfb(start, resfn, jacfn=NULL, trace=FALSE, lower=-Inf, upper=Inf, maskidx=NULL, control, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resfn}{
     A function that evaluates the residual vector for computing the elements of
     the sum of squares function at the set of parameters \code{start}. 
}
  \item{jacfn}{
     A function that evaluates the Jacobian of the sum of squares function, that
     is, the matrix of partial derivatives of the residuals with respect to each
     of the parameters. If NULL (default), uses an approximation.
     ?? put in character form as in optimx??
}
  \item{start}{
     A named parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     \code{nls()} takes a list, and that is permitted here also.
}
  \item{trace}{
     Logical TRUE if we want intermediate progress to be reported. Default is FALSE.
} 
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default -Inf.
} 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default Inf.
} 
  \item{maskidx}{
     Vector if indices of the parameters to be masked. These parameters will NOT be altered
     by the algorithm. Note that the mechanism here is different from that in \code{nlxb}
     which uses the names of the parameters.

} 
  \item{control}{ 
     A list of controls for the algorithm. These are:
%%     ~~Describe \code{control} here~~
}
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same. It is not clear how the dot variables should
     be used, since data should be in 'data'.
}
}
\details{
    \code{nlfb} attempts to solve the nonlinear sum of squares problem by using
    a variant of Marquardt's approach to stabilizing the Gauss-Newton method using
    the Levenberg-Marquardt adjustment. This is explained in Nash (1979 or 1990) in
    the sections that discuss Algorithm 23. (?? do we want a vignette. Yes, because
    folk don't have access to book easily, but finding time.)

    In this code, we solve the (adjusted) Marquardt equations by use of the 
    \code{qr.solve()}. Rather than forming the J'J + lambda*D matrix, we augment
    the J matrix with extra rows and the y vector with null elements. 
}
\value{
  A list of the following items
  \item{coeffs}{A named vector giving the parameter values at the supposed solution.}
  \item{ssquares}{The sum of squared residuals at this set of parameters.}
  \item{resid}{The residual vector at the returned parameters.}
  \item{jacobian}{The jacobian matrix (partial derivatives of residuals w.r.t. the parameters)
           at the returned parameters.}
  \item{feval}{The number of residual evaluations (sum of squares computations) used.}
  \item{jeval}{The number of Jacobian evaluations used.}
}
\references{
%% ~put references to the literature/web site here ~
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!

}
\author{
John C Nash <nashjc@uottawa.ca>
}
\note{
   Special notes, if any, will appear here.
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}

\examples{
cat("See examples in nlmrt-package.Rd\n")

}
\keyword{ nonlinear least squares }

