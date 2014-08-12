\name{modeldoc}
\alias{modeldoc}
\title{
Output a description of the model and the data that was used at the time
of its creation to the console and optionally to a file. The purpose of 
this function is to provide a record of the details underlying the function
\code{fun} as well as to aid users wishing to create objective functions for
optimization tools.
}
\description{
These functions create functions to evaluate residuals or sums of squares at 
particular parameter locations.
}
\usage{
modeldoc(fun, savefile=NULL)
}
\arguments{
  \item{fun}{
A function produced by one of \code{model2rjfun} or \code{model2ssgrfun}.
}

  \item{savefile}{A character string giving the filename of a file
to record the terminal output. The filename may include a path relative
to the working directory.
}

}
\details{
\code{modeldoc} displays the contents of the environment associated with
\code{fun}. 

**WARNING** \code{modeldoc} currently lists variables from the \code{"..."}
arguments but does NOT seem to be able to list the values of these variables
that are current at the time \code{modeldoc} is called.
}

\value{
The function returns 0 unless it stops because \code{savefile} is 
neither NULL nor a character string.
}

\author{
John Nash and Duncan Murdoch
}

\seealso{
\code{\link{nls}}
}
\examples{

## require(nls14)
traceval  <-  TRUE  # traceval set TRUE to debug or give full history
# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing
# A simple starting vector -- must have named parameters for nls14xb, nls, wrap14nls.
start1  <-  c(b1=1, b2=1, b3=1)
weeddata2  <-  data.frame(y=1.5*ydat, tt=tdat)
escal <-  y ~ ms*b1/(1+b2*exp(-b3*tt))
ms <- 1
weeddata3<-weeddata2
weeddata3$y<-0
grs3<-model2rjfunx(escal, start1, data=weeddata3)
res3<-grs3(start1)
res3
modeldoc(grs3)

# Now a different value of ms
ms<-2
grs3b<-model2rjfunx(escal, start1, data=weeddata3)
res3b<-grs3b(start1)
res3b
modeldoc(grs3b)
# modeldoc(grs3b, savefile="grs3save.txt") ## to save the output

}
\keyword{ nonlinear }
