\name{nlsmnqb}
\alias{nlsmnqb}
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
   nlsmnqb(formula, start, trace=FALSE, data, lower=-Inf, upper=Inf, masked=NULL, control, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
     This is a modeling formula of the form (as in \code{nls})
           lhsvar ~ rhsexpression
     for example,
           y ~ b1/(1+b2*exp(-b3*T))
     You may also give this as a string.
}
  \item{start}{
     A named parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     \code{nls()} takes a list, and that is permitted here also.
}
  \item{trace}{
     Logical TRUE if we want intermediate progress to be reported. Default is FALSE.
} 
  \item{data}{
     A data frame containing the data of the variables in the formula. This data may,
     however, be supplied directly in the parent frame. 
} 
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default -Inf.
} 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default Inf.
} 
  \item{masked}{
     Character vector of quoted parameter names. These parameters will NOT be altered
     by the algorithm.
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
    \code{nlsmnqb} attempts to solve the nonlinear sum of squares problem by using
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

ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat<-1:length(ydat) # for testing
# WARNING -- using T can get confusion with TRUE
start1<-c(b1=1, b2=1, b3=.1)
eunsc<- y ~ b1/(1+b2*exp(-b3*t))
weeddata<-data.frame(y=ydat, t=tdat)

anlsmnqb1<-try(nlsmn0b(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b0=500, b2=100, b3=5), trace=TRUE, data=weeddata))
print(anlsmnqb1)

tmp<-readline("next")

anlsb1<-try(nls(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b0=500, b2=100, b3=5), trace=TRUE, data=weeddata, algorithm='port'))
print(anlsb1)

tmp<-readline("next")

anlsmnqb2<-try(nlsmn0b(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b0=500, b2=100, b3=.25), trace=TRUE, data=weeddata))
print(anlsmnqb2)

tmp<-readline("next")

anlsb2<-try(nls(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b0=500, b2=100, b3=.25), trace=TRUE, data=weeddata, algorithm='port'))
print(anlsb2)

tmp<-readline("next")

cat("TEST MASKING OF b2\n")

anlsmnqm<-try(nlsmn0b(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), upper=c(b0=500, b2=100, b3=5), masked=c("b2"), trace=TRUE, data=weeddata))
print(anlsmnqm)

tmp<-readline("next")

anlsmnqmx<-try(nlsmn0b("y ~ b1/(1+exp(-b3*t))", start=c(b1=1, b3=.1), lower=c(b1=0, b3=0), upper=c(b0=500, b3=.25), trace=TRUE, data=weeddata))
print(anlsmnqmx)

tmp<-readline("next")


ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat<-1:length(ydat) # for testing
wedat<-data.frame(t=tdat, y=ydat)
# WARNING -- using T can get confusion with TRUE
start1<-c(b1=1, b2=1, b3=1)
eunsc<- y ~ b1/(1+b2*exp(-b3*t))

cat("UNCONSTRAINED\n")
an1q<-try(nlsmnqb(eunsc, start=start1, trace=TRUE, data=wedat))
print(an1q)
tmp<-readline("next")


cat("MASKED\n")
an1qm3<-try(nlsmnqb(eunsc, start=start1, trace=TRUE, data=wedat, masked=c("b3")))
print(an1qm3)
tmp<-readline("next")

# Note that the parameters are put in out of order to test code.
an1qm123<-try(nlsmnqb(eunsc, start=start1, trace=TRUE, data=wedat, masked=c("b2","b1","b3")))
print(an1qm123)
tmp<-readline("next")


cat("BOUNDS")
start2<-c(b1=100, b2=10, b3=0.1)
an1qb1<-try(nlsmnqb(eunsc, start=start2, trace=TRUE, data=wedat, lower=c(0,0,0), upper=c(200, 60, .3)))
print(an1qb1)
tmp<-readline("next")


cat("BOUNDS and MASK")
an1qbm2<-try(nlsmnqb(eunsc, start=start2, trace=TRUE, data=wedat, lower=c(0,0,0), upper=c(200, 60, .3), masked=c("b2")))
print(an1qbm2)
tmp<-readline("next")



}
\keyword{ nonlinear least squares }

