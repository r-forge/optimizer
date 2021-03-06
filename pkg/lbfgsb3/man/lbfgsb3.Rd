\name{lbfgsb3}
\alias{lbfgsb3}
\title{Interfacing wrapper for the Nocedal - Morales LBFGSB3 (Fortran) 
      limited memory BFGS solver.
}
\description{
   This package is an 
}
\usage{
   lbfgsb3(prm, fn, gr=NULL, lower=-Inf, upper=Inf, 
         control=list(), ...)
}
\arguments{
  \item{prm}{
     A parameter vector which gives the initial guesses to the parameters
     that will minimize \code{fn}. This can be named, for example, we could use
        prm=c(b1=1, b2=2.345, b3=0.123)
  }
  \item{fn}{
     A function that evaluates the objective function to be minimized.
  }
  \item{gr}{
     If present, a function that evaluates the gradient vector for 
     the objective function at the given parameterscomputing the elements of
     the sum of squares function at the set of parameters \code{start}. 
  }
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default -Inf.
  } 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default Inf.
  } 
  \item{control}{ An optional list of control settings. See below in details.}
  \item{\dots}{Any data needed for computation of the objective function and gradient.}
     
}
\details{

  See the notes below for a general appreciation of this package.

  The list of controls for the algorithm contains:
  \describe{
      \item{\code{trace}}{an integer which if 0 (default) causes no intermediate
         output in the R section of the code, otherwise some diagnostic information.
         The larger the value of \code{trace} the greater the amount of information.
         Because output from Fortran does not play well with R, much of the 
         diagnostic output of the Fortran has been suppressed.}
      \item{\code{maxit}}{Maximum number of function/gradient evaluations. }
  }


  The output items include the following (not fully edited at 150121).

  info <- list(task = task, itask = itask, lsave = lsave, 
                icsave = icsave, dsave = dsave, isave = isave)

      icsave is a working integer

      lsave is a logical working array of dimension 4.
        On exit with 'task' = NEW_X, the following information is 
                                                              available:
          If lsave(1) = .true.  then  the initial X has been replaced by
                                      its projection in the feasible set;
          If lsave(2) = .true.  then  the problem is constrained;
          If lsave(3) = .true.  then  each variable has upper and lower
                                      bounds;

      isave is an integer working array of dimension 44.
        On exit with 'task' = NEW_X, the following information is 
                                                              available:
          isave(22) = the total number of intervals explored in the 
                          search of Cauchy points;
          isave(26) = the total number of skipped BFGS updates before 
                          the current iteration;
          isave(30) = the number of current iteration;
          isave(31) = the total number of BFGS updates prior the current
                          iteration;
          isave(33) = the number of intervals explored in the search of
                          Cauchy point in the current iteration;
          isave(34) = the total number of function and gradient 
                          evaluations;
          isave(36) = the number of function value or gradient
                                   evaluations in the current iteration;
          if isave(37) = 0  then the subspace argmin is within the box;
          if isave(37) = 1  then the subspace argmin is beyond the box;
          isave(38) = the number of free variables in the current
                          iteration;
          isave(39) = the number of active constraints in the current
                          iteration;
          n + 1 - isave(40) = the number of variables leaving the set of
                            active constraints in the current iteration;
          isave(41) = the number of variables entering the set of active
                          constraints in the current iteration.

      dsave is a double precision working array of dimension 29.
        On exit with 'task' = NEW_X, the following information is
                                                              available:
          dsave(1) = current 'theta' in the BFGS matrix;
          dsave(2) = f(x) in the previous iteration;
          dsave(3) = factr*epsmch;
          dsave(4) = 2-norm of the line search direction vector;
          dsave(5) = the machine precision epsmch generated by the code;
          dsave(7) = the accumulated time spent on searching for
                                                          Cauchy points;
          dsave(8) = the accumulated time spent on
                                                  subspace minimization;
          dsave(9) = the accumulated time spent on line search;
          dsave(11) = the slope of the line search function at
                                   the current point of line search;
          dsave(12) = the maximum relative step length imposed in
                                                            line search;
          dsave(13) = the infinity norm of the projected gradient;
          dsave(14) = the relative step length in the line search;
          dsave(15) = the slope of the line search function at
                                  the starting point of the line search;
          dsave(16) = the square of the 2-norm of the line search
                                                       direction vector.

}
\value{
  A list of the following items
  \item{prm}{A vector giving the parameter values at the supposed solution.}
  \item{f}{The value of the objective function at this set of parameters.}
  \item{g}{An estimate of the gradient of the objective at the solution.}
  \item{info}{A structure containing information on the disposition of the
       computation on return. See Details.}
}
\references{
     Morales, J. L.; Nocedal, J. (2011). "Remark on 'algorithm 778: L-BFGS-B: 
           Fortran subroutines for large-scale bound constrained optimization' ". 
           ACM Transactions on Mathematical Software 38: 1. 

     Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C. (1995). "A Limited Memory Algorithm 
           for Bound Constrained Optimization". SIAM J. Sci. Comput. 16 (5): 1190-1208. 
           
     Zhu, C.; Byrd, Richard H.; Lu, Peihuang; Nocedal, Jorge (1997). "L-BFGS-B: 
           Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained 
           optimization". ACM Transactions on Mathematical Software 23 (4): 550-560.
}
\author{
     John C Nash <nashjc@uottawa.ca> (of the wrapper and edits to Fortran code to allow R output)
     Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales (original Fortran packages)     
}
\note{
   This package is a wrapper to the Fortran code released by Nocedal and Morales.
   This poses several difficulties for an R package. While the \code{.Fortran()}
   tool exists for the interfacing, we must be very careful to align the arguments
   with those of the Fortran subroutine, especially in type and storage. 
   
   A more annoying task for interfacing the Fortran code is that Fortran WRITE or
   PRINT statements must all be replaced with calls to special R-friendly output 
   routines. Unfortunately, the Fortran is full of output statements. Worse, we may
   wish to be able to suppress such output, and there are thus many modifications
   to be made. This means that an update of the original code cannot be simply
   plugged into the R package \code{src} directory.
   
   Finally, and likely because L-BFGS-B has a long history, the Fortran code is far
   from well-structured. For example, the number of function and gradient evaluations
   used is returned as the 34'th element of an integer vector. There does not appear
   to be an easy way to stop the program after some maximum number of such evaluations
   have been performed.

   On the other hand, the version of L-BFGS-B in \code{optim()} is a \code{C} translation
   of a now-lost Fortran code. It does not implement the improvements Nocedal and 
   Morales published in 2011. Hence, despite its deficiencies, this wrapper has been
   prepared.

}

\seealso{
    Packages \code{\link{optim}} and \code{optimx}.
}

\examples{
cat("Examples are to be added\n")
cat("But see the tests\n")

}
\keyword{ nonlinear parameter optimization }

