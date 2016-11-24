\name{runoptprob}
\alias{runoptprob}
\title{General optimization test runner}
\concept{minimization}
\concept{maximization}
\concept{computational testing}
\description{
  A tool to run various pre-defined problems with many different 
  optimization or nonlinear least squares methods. 
}
\usage{
runoptprob <- function(pfilename, probclass=NULL, minmeth=NULL, submeth=NULL, istart="1", 
                 runopts=list(), control=list(), ...) 
}
\arguments{
 \item{pfilename}{a character string giving the filename of a particular problem located in the 
       inst/extdata/ directory of the package.}  
 \item{probclass}{If NULL (default) then attempt all available classes, else one of "uncopt",
       "boundopt", "formula", or "sumsquares".}
 \item{minmeth}{One of a selected set of minimization tool classes. The list includes??.}
 \item{submeth}{A particular member of the class of methods in minmeth.}   
 \item{istart}{The character index of a pre-defined starting set of parameters for
     the specified problem file. This is a character index to allow for special case
     e.g., "10:pi" would be an index for the starting set of parameters of length 10 
    each of value pi, which is a reasonable choice for an extended Rosenbrock function.}
 \item{runopts}{Any arguments for running a particular minimization method needed.}
 \item{control}{Any controls that should be supplied to the minimization method.}
 \item{\dots}{Any other (named) arguments}
}
\details{
  To be filled in!??
}
\value{

  A list with components:
  \describe{
   \item{par}{The best set of parameters found.}
   \item{pfilename}{The name of the original problem file (root thereof??).}
   \item{nstart}{ The index of the starting parameter set.}
   \item{minmeth}{ The class of the minimization method. }
   \item{submeth}{ The particular member of the class of minimizers.}
  }
}
\keyword{nonlinear}
\keyword{optimize}
