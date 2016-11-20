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
runoptprob <- function(pfilename, minmeth=NULL, submeth=NULL, nstart=0, 
                 runopts=list(), control=list(), ...) 
}
\arguments{
 \item{pfilename}{a character string giving the filename of a particular problem located in the inst/extdata/ directory of the package.}  
 \item{minmeth}{One of a selected set of minimization tool classes. The list includes??.}
 \item{submeth}{A particular member of the class of methods in minmeth.}   
 \item{nstart}{The index of the pre-defined starting set of parameters in the specified problem file, or 0 for all predefined starting vectors, ?? or (?? could have user defined if negative??)}
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
