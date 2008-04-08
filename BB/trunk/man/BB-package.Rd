\name{BB-package}
\alias{BB-package}
\alias{BB.Intro}
\docType{package}

\title{Solving and Optimizing Large-Scale Nonlinear Problems}

\description{Non-monotone Barzilai-Borwein spectral methods for the solution and optimization of large-scale nonlinear systems.  }

\details{
\tabular{ll}{
Package: \tab BB\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 3.\cr
}
The main functions are:
\preformatted{

dfsane  function for solving large systems of nonlinear equations using a 
        derivative-free spectral approach 

sane	function for solving large systems of nonlinear equations using 
        spectral approach

spg	function for spectral projected gradient method for large-scale 
        optimization with simple constraints
}

Maintainer:Ravi Varadhan <rvaradhan@jhmi.edu> and 
   Paul Gilbert <pgilbert@bank-banque-canada.ca> 
}
\author{Ravi Varadhan}
\references{

   J Barzilai, and JM Borwein  (1988),  Two-point step size gradient methods, \emph{IMA J Numerical Analysis}, 8, 141-148. 

    Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected 
    gradient methods on convex sets, \emph{SIAM J Optimization}, 10, 1196-1211.
    
    Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for 
    convex-constrained optimization, \emph{ACM Transactions on Mathematical Software}.
   
   L Grippo, F Lampariello, and S Lucidi (1986),  A nonmonotone line search technique for Newton's method,  \emph{SIAM J on Numerical Analysis}, 23, 707-716. 

   W LaCruz, and M Raydan (2003),  Nonmonotone spectral methods for large-scale nonlinear systems, \emph{Optimization Methods and Software}, 18, 583-599. 

   W LaCruz, JM Martinez, and M Raydan (2006),  Spectral residual mathod without gradient information for solving large-scale nonlinear systems of equations, \emph{Mathematics of Computation}, 75, 1429-1448. 

   M Raydan (1997),  Barzilai-Borwein gradient method for large-scale unconstrained minimization problem, \emph{SIAM J of Optimization}, 7, 26-33.

   R Varadhan and C Roland (2008),  Simple and globally-convergent methods for accelerating the convergence of any EM 
   algorithm,  \emph{Scandinavian J Statistics}. 

   R Varadhan and P Gilbert (2008),  BB: A package of Barzilai-Borwein spectral methods for solving and optimizing large-scale nonlinear systems, Unpublished.

}
\keyword{package}

