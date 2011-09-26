\name{clg}
\alias{clg}
\alias{PROBj2k}
\alias{GrPROBj2k}
\alias{GrGrPROBj2k}
\title{Probability to lapse}
\description{
		Cumulative logistic function and its successive derivatives
	}
\usage{
PROBj2k(x, j, k, param_j)
GrPROBj2k(x, j, k, ideriv, param_j)
GrGrPROBj2k(x, j, k, ideriv, mderiv, param_j)
}
\arguments{
  \item{x}{vector \code{x}.}
  \item{j}{integer \code{j} for player 'from'.}
  \item{k}{integer \code{k} for player 'to'.}
  \item{ideriv}{integer \code{ideriv} for first derivative index.}
  \item{mderiv}{integer \code{mderiv} for second derivative index.}
  \item{param_j}{vector \code{param_j}.}    
}
\value{
A numeric.
}
\author{
 Christophe Dutang
}
\keyword{math}


