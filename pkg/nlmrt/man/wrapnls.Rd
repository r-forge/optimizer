\name{wrapnls}
\alias{wrapnls}
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
   wrapnls(formula, start, trace=FALSE, data, control, ...)
}
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
}
  \item{trace}{
     Logical TRUE if we want intermediate progress to be reported. Default is FALSE.
} 
  \item{data}{
     A data frame containing the data of the variables in the formula. This data may,
     however, be supplied directly in the parent frame. 
} 
  \item{control}{ 
     A list of controls for the algorithm. These are:
%%     ~~Describe \code{control} here~~
}
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{wrapnls} first attempts to solve the nonlinear sum of squares problem by using
    \code{nlsmnq}, then takes the parameters from that method to call \code{nls}.
}
\value{
  An object of type nls.}
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
start1<-c(b1=100, b2=10, b3=1)
eunsc<- y ~ b1/(1+b2*exp(-b3*t))
wdata<-data.frame(t=tdat, y=ydat)

an1w<-try(wrapnls(eunsc, start=start1, trace=TRUE, data=wdata))
print(an1w)

## Put in example using data frame!
}
\keyword{ nonlinear least squares }

