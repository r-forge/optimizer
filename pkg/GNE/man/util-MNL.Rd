\name{MNL}
\alias{MNL}
\alias{PROBj2kRatio}
\alias{GrPROBj2kRatio}
\alias{GrGrPROBj2kRatio}
\alias{PROBj2kDiff}
\alias{GrPROBj2kDiff}
\alias{GrGrPROBj2kDiff}
\title{Probability to lapse from j to k.}
\description{
		Multinomial logit function and its successive derivatives
	}
\usage{
PROBj2kRatio(x, j, k, param_j)
GrPROBj2kRatio(x, j, k, ideriv, param_j)
GrGrPROBj2kRatio(x, j, k, ideriv, mderiv, param_j)
PROBj2kDiff(x, j, k, param_j)
GrPROBj2kDiff(x, j, k, ideriv, param_j)
GrGrPROBj2kDiff(x, j, k, ideriv, mderiv, param_j)
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


