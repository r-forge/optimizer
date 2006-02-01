\name{00.GPArotation.Intro}
\alias{00.GPArotation.Intro}
\alias{GPArotation.Intro}
\docType{package}

\title{GPA Rotation for Factor Analysis}

\description{GPArotation implements Gradient Projection Algorithms and several
rotation objective functions for factor analysis.
}

\details{
\tabular{ll}{
Package: \tab GPArotation\cr
Depends: \tab R (>= 2.0.0)\cr
License: \tab GPL Version 2.\cr
URL: \tab http://www.stat.ucla.edu/research or\cr
     \tab http://www.stat.ucla.edu/research/gpa\cr}

The main optimization functions are
\code{\link{GPForth}} and \code{\link{GPFoblq}}. Rotation objectives include
\code{\link{oblimin}} and many others.
}

\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert.

Code is modified from original source \file{splusfunctions.net} available at 
  \url{http://www.stat.ucla.edu/research/gpa}.

Maintainer: Paul Gilbert <pgilbert@bank-banque-canada.ca>  
}
\references{
   The software reference is

   Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

   Theory of gradient projection algorithms may be found in:

   Jennrich, R.I. (2001). A simple general procedure for orthogonal rotation.
       \emph{Psychometrika}, \bold{66}, 289--306.

   Jennrich, R.I. (2002). A simple general method for oblique rotation. 
       \emph{Psychometrika}, \bold{67}, 7--19.
}
\keyword{ package }

\seealso{
   \code{\link{GPForth}},
   \code{\link{GPFoblq}},
   \code{\link[stats]{factanal}} 
}
