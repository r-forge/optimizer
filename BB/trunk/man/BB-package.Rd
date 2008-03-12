\name{BB-package}
\alias{BB-package}
\alias{BB.Intro}
\docType{package}

\title{Barzilai-Borwein Optimization}

\description{Large-scale optimization and root-finding.}

\details{
\tabular{ll}{
Package: \tab BB\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 3.\cr
}
The main functions are
\preformatted{
spg	  functions for spectral projected gradient method for
  	  box-constrained minimization

sane	  functions for spectral-gradient based solution of nonlinear systems

dfsane	  functions for spectral-gradient based solution of
  	  nonlinear systems (derivative-free sane algorithm)

nlsolve	  function that solves a nonlinear system by converting
  	  it into a minimization problem and using optim's L-BFGS-B algorithm
}

Maintainer:Ravi Varadhan <rvaradhan@jhmi.edu> and 
   Paul Gilbert <pgilbert@bank-banque-canada.ca> 
}
\author{Ravi Varadhan}
\references{
    Non-monotone spectral projected-gradient method for minimization
    Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected 
    gradient methods on convex sets, \emph{SIAM J Optimization}, 10, 1196-1211.
    
    Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for 
    convex-constrained optimization, 
    \emph{ACM Transactions on Mathematical Software}.

    LaCruz and Raydan M (2003): 

    Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
}
\keyword{package}

