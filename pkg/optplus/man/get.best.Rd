\name{get.best}
\alias{get.best}
\encoding{UTF-8}
\title{Return best optimx result in the form that optim() returns}
\concept{minimization}
\concept{maximization}
\description{
  The optimx general-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing \code{optim()} function, creates
  a data frame of results. \code{get.best} selects the best of these (at time of
  writing maximization problems may not work) and returns this answer in the form
  that matches the return from \code{optim}.
}
\usage{
   get.best(optimx.obj, maximize=FALSE) 
}
\arguments{
 \item{optimx.obj}{object returned by `optimx'}
 \item{maximize}{TRUE if the objective function is to be maximized.}
}
\details{
  None
}
\value{
  \code{get.best} returns a list
  \item{par}{The vector of parameters at the supposed optimum}
  \item{value}{The value of the objective at the supposed optimum}
  \item{counts}{A vector of length 2 with the number of function and 
       gradient calls respectively.}
  \item{convergence}{A code indicating the status of the result. See relevant optimizer.}
  \item{message}{A message, if any, returned by the optimizer. To allow the user to know
     which method was used to obtain this result, the method is appended to this message.}
  \item{hessian}{If available, the (approximate) Hessian matrix at the returned parameters.}
}
\keyword{nonlinear}
\keyword{optimize}

