##' @importFrom Rcpp evalCpp
## below added
##' @importFrom Rcpp sourceCpp
##' @importFrom methods is
##' @useDynLib lbfgsb3x, .registration=TRUE
"lbfgsb3"

##' Interfacing wrapper for the Nocedal - Morales LBFGSB3 (Fortran) 
##'       limited memory BFGS solver.
##' 
##' @param prm
##'      A parameter vector which gives the initial guesses to the parameters
##'      that will minimize \code{fn}. This can be named, for example, we could use
##'         prm=c(b1=1, b2=2.345, b3=0.123)
##' @param fn
##'      A function that evaluates the objective function to be minimized.
##' @param  gr
##'      If present, a function that evaluates the gradient vector for 
##'      the objective function at the given parameterscomputing the elements of
##'      the sum of squares function at the set of parameters \code{start}. 
##' @param  lower
##'      Lower bounds on the parameters. If a single number, this will be applied to all
##'      parameters. Default -Inf.
##' @param  upper
##'      Upper bounds on the parameters. If a single number, this will be applied to all
##'      parameters. Default Inf.
##' @param  control An optional list of control settings. See below in details.
##' @param  \dots Any data needed for computation of the objective function and gradient.
##' 
##' 
##' @details
##' 
##' See the notes below for a general appreciation of this package.
##'
##' The control list can contain:
##' \itemize{
##'    \item{\code{trace}}{an integer which if 0 (default) causes no intermediate
##'       output in the R section of the code, otherwise some diagnostic information.
##'       The larger the value of \code{trace} the greater the amount of information.
##'       In the Fortran code, see the variable \code{iprint}, which is equal to \code{trace} for
##'       positive values of \code{trace}, and set to -1 for \code{trace=0}. (\code{iprint=0} is not 
##'       possible via this R wrapper.)}
##'    \item{\code{maxit}}{Maximum number of function/gradient evaluations. }
##'    \item{info} {a boolean to indicate if more optimization information is captured
##'     and output in a \code{info} list}
##' 
##' 
##'   The output items include the following (not fully edited at 150121).
##' 
##'   info <- list(task = task, itask = itask, lsave = lsave, 
##'                 icsave = icsave, dsave = dsave, isave = isave)
##' 
##'       icsave is a working integer
##' 
##'       lsave is a logical working array of dimension 4.
##'         On exit with 'task' = NEW_X, the following information is 
##'                                                               available:
##'           If lsave(1) = .true.  then  the initial X has been replaced by
##'                                       its projection in the feasible set;
##'           If lsave(2) = .true.  then  the problem is constrained;
##'           If lsave(3) = .true.  then  each variable has upper and lower
##'                                       bounds;
##' 
##'       isave is an integer working array of dimension 44.
##'         On exit with 'task' = NEW_X, the following information is 
##'                                                               available:
##'           isave(22) = the total number of intervals explored in the 
##'                           search of Cauchy points;
##'           isave(26) = the total number of skipped BFGS updates before 
##'                           the current iteration;
##'           isave(30) = the number of current iteration;
##'           isave(31) = the total number of BFGS updates prior the current
##'                           iteration;
##'           isave(33) = the number of intervals explored in the search of
##'                           Cauchy point in the current iteration;
##'           isave(34) = the total number of function and gradient 
##'                           evaluations;
##'           isave(36) = the number of function value or gradient
##'                                    evaluations in the current iteration;
##'           if isave(37) = 0  then the subspace argmin is within the box;
##'           if isave(37) = 1  then the subspace argmin is beyond the box;
##'           isave(38) = the number of free variables in the current
##'                           iteration;
##'           isave(39) = the number of active constraints in the current
##'                           iteration;
##'           n + 1 - isave(40) = the number of variables leaving the set of
##'                             active constraints in the current iteration;
##'           isave(41) = the number of variables entering the set of active
##'                           constraints in the current iteration.
##' 
##'       dsave is a double precision working array of dimension 29.
##'         On exit with 'task' = NEW_X, the following information is
##'                                                               available:
##'           dsave(1) = current 'theta' in the BFGS matrix;
##'           dsave(2) = f(x) in the previous iteration;
##'           dsave(3) = factr*epsmch;
##'           dsave(4) = 2-norm of the line search direction vector;
##'           dsave(5) = the machine precision epsmch generated by the code;
##'           dsave(7) = the accumulated time spent on searching for
##'                                                           Cauchy points;
##'           dsave(8) = the accumulated time spent on
##'                                                   subspace minimization;
##'           dsave(9) = the accumulated time spent on line search;
##'           dsave(11) = the slope of the line search function at
##'                                    the current point of line search;
##'           dsave(12) = the maximum relative step length imposed in
##'                                                             line search;
##'           dsave(13) = the infinity norm of the projected gradient;
##'           dsave(14) = the relative step length in the line search;
##'           dsave(15) = the slope of the line search function at
##'                                   the starting point of the line search;
##'           dsave(16) = the square of the 2-norm of the line search
##'                                                        direction vector.
##' 
##' }
##'
##' @return
##'   A list of the following items
##' \itemize{
##' \item{prm} The best set of parameters found.
##' \item{f} The value of fn corresponding to prm.
##' \item{g}{An estimate of the gradient of the objective at the solution.}
##'   \item{info}{A structure containing information on the disposition of the
##'        computation on return. See Details.}
##' }
##' @seealso Packages \code{\link{optim}} and \code{optimx}.
##' @keywords nonlinear parameter optimization
##' @author Matthew Fidler (move to C and add more options for adjustments),
##'     John C Nash <nashjc@uottawa.ca> (of the wrapper and edits to Fortran code to allow R output)
##'     Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales (original Fortran packages)
##'
##' @references
##'     Morales, J. L.; Nocedal, J. (2011). "Remark on 'algorithm 778: L-BFGS-B:
##'           Fortran subroutines for large-scale bound constrained optimization' ".
##'           ACM Transactions on Mathematical Software 38: 1.
##'
##'     Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C. (1995). "A Limited Memory Algorithm
##'           for Bound Constrained Optimization". SIAM J. Sci. Comput. 16 (5): 1190-1208.
##'
##'     Zhu, C.; Byrd, Richard H.; Lu, Peihuang; Nocedal, Jorge (1997). "L-BFGS-B:
##'           Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained
##'           optimization". ACM Transactions on Mathematical Software 23 (4): 550-560.
##'
##' @note
##'   This package is a wrapper to the Fortran code released by Nocedal and Morales.
##'   This poses several difficulties for an R package. While the \code{.Fortran()}
##'   tool exists for the interfacing, we must be very careful to align the arguments
##'   with those of the Fortran subroutine, especially in type and storage.
##'
##'   A more annoying task for interfacing the Fortran code is that Fortran WRITE or
##'   PRINT statements must all be replaced with calls to special R-friendly output
##'   routines. Unfortunately, the Fortran is full of output statements. Worse, we may
##'   wish to be able to suppress such output, and there are thus many modifications
##'   to be made. This means that an update of the original code cannot be simply
##'   plugged into the R package \code{src} directory.
##'
##'   Finally, and likely because L-BFGS-B has a long history, the Fortran code is far
##'   from well-structured. For example, the number of function and gradient evaluations
##'   used is returned as the 34'th element of an integer vector. There does not appear
##'   to be an easy way to stop the program after some maximum number of such evaluations
##'   have been performed.
##'
##'   On the other hand, the version of L-BFGS-B in \code{optim()} is a \code{C} translation
##'   of a now-lost Fortran code. It does not implement the improvements Nocedal and
##'   Morales published in 2011. Hence, despite its deficiencies, this wrapper has been
##'   prepared.
##'
##' In addition to the above reasons for the original lbfgsb3 package,
##' this additional package allows C calling of L-BFGS-B 3.0 by a
##' program as well as adjustments to the tolerances that were not
##' present in the original CRAN package.  Also adjustments were made
##' to have outputs conform with R's optim routine.
##' @examples
##' # Rosenbrock's banana function
##' n=3; p=100
##'
##' fr = function(x)
##' {
##'     f=1.0
##'     for(i in 2:n) {
##'         f=f+p*(x[i]-x[i-1]**2)**2+(1.0-x[i])**2
##'     }
##'     f
##' }
##'
##' grr = function(x)
##' {
##'     g = double(n)
##'     g[1]=-4.0*p*(x[2]-x[1]**2)*x[1]
##'     if(n>2) {
##'         for(i in 2:(n-1)) {
##'             g[i]=2.0*p*(x[i]-x[i-1]**2)-4.0*p*(x[i+1]-x[i]**2)*x[i]-2.0*(1.0-x[i])
##'         }
##'     }
##'     g[n]=2.0*p*(x[n]-x[n-1]**2)-2.0*(1.0-x[n])
##'     g
##' }
##' x = c(a=1.02, b=1.02, c=1.02)
##' (op1 <- lbfgsb3c(x,fr, grr, x))
##' @export
lbfgsb3 <- function(prm, fn, gr=NULL, lower = -Inf, upper = Inf,
         control=list(), ...){
# ?? need to add controls
# if (is.null(gr)) require(numDeriv) # eventually change to "grfwd" etc.
# interface to Fortran Lbfgsb.3.0
## 150121 There is currently no limit on function or gradient evaluations

    tasklist <- c('NEW_X', 'START', 'STOP', 'FG',  # 1-4
       'ABNORMAL_TERMINATION_IN_LNSRCH', 'CONVERGENCE', #5-6
       'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL',#7
       'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH',#8
       'ERROR: FTOL .LT. ZERO', #9
       'ERROR: GTOL .LT. ZERO' ,#10
       'ERROR: INITIAL G .GE. ZERO', #11
       'ERROR: INVALID NBD', # 12
       'ERROR: N .LE. 0', # 13
       'ERROR: NO FEASIBLE SOLUTION', # 14
       'ERROR: STP .GT. STPMAX', # 15
       'ERROR: STP .LT. STPMIN', # 16
       'ERROR: STPMAX .LT. STPMIN', # 17
       'ERROR: STPMIN .LT. ZERO', # 18
       'ERROR: XTOL .LT. ZERO', # 19
       'FG_LNSRCH', # 20
       'FG_START', # 21
       'RESTART_FROM_LNSRCH', # 22
       'WARNING: ROUNDING ERRORS PREVENT PROGRESS', # 23
       'WARNING: STP .eq. STPMAX', # 24
       'WARNING: STP .eq. STPMIN', # 25
       'WARNING: XTOL TEST SATISFIED')# 26
# CONV in 6, 7, 8; ERROR in 9-19; WARN in 23-26
 
# if (!is.loaded("lbfgsb3.so")) dyn.load("lbfgsb3.so") # get the routines attached

    factr <- 1.0e+7
    pgtol <- 1.0e-5
    nmax <- 1024L
    mmax <- 17L

    if (length(prm) > nmax) stop("The number of parameters cannot exceed 1024")
    n <- as.integer(length(prm))
    m <- 5L # default 

# control defaults -- idea from spg
ctrl <- list(trace = 0, maxit = 100*n) ## ??  iprint = 0L)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control

if (ctrl$trace > 2) print(control)


# Here expand control list, but for moment leave alone
      iprint <- as.integer(ctrl$trace)
      if (iprint == 0) { iprint <- as.integer(-1L) } # Note this excludes iprint=0 case
      # from working in the Fortran code, but because of R restrictions on Fortran
      # output, some intermediate output is suppressed.

## Define the storage
nbd <- rep(2L, n) # start by defining them "on" -- adjust below
nwa<-2*mmax*nmax + 5*nmax + 11*mmax*mmax + 8*mmax
wa<-rep(0, nwa)
dsave<-rep(0,29)
lsave<-rep(TRUE,4)
isave<-rep(0L,44)
iwa<-rep(0L, 3*nmax)
csave<-"" # note char strings are 255 automatically


if (length(lower) == 1) lower <- rep(lower, n)
if (length(upper) == 1) upper <- rep(upper, n)

bigval <- .Machine$double.xmax/10.

for (i in 1:n) {
   if (is.finite(lower[i])) {
        if (is.finite(upper[i])) nbd[i] <- 2
        else {
           nbd[i] <- 1
           upper[i] <- bigval # to avoid call issue
             }
   } else { if (is.finite(upper[i])) {
              nbd[i] <- 3
              lower[i] <- -bigval
            } else {
              nbd[i] <- 0 
              upper[i] <- bigval
              lower[i] <- -bigval
                 }
   }
}
## cat("nbd:")
## print(nbd)


##     We start the iteration by initializing task.
## 
      itask <- 2L # START
      task <- tasklist[itask]
      f <- .Machine$double.xmax / 100
      g <- rep(f, n)
      nfg <- 0 # initialize counter

##        ------- the beginning of the loop ----------
icsave <- 0 # to make sure defined

## 111  continue ##  top of loop
repeat {
##     This is the call to the L-BFGS-B code.
      if (ctrl$trace >= 2) {
       cat("Before call, f=",f,"  task number ",itask," ")
       print(task)
      }
      result <- .Fortran('lbfgsb3', n = as.integer(n),m = as.integer(m),
                   x = as.double(prm), l = as.double(lower), u = as.double(upper),
                   nbd = as.integer(nbd), f = as.double(f), g = as.double(g),
                   factr = as.double(factr), pgtol = as.double(pgtol),
                   wa = as.double(wa), iwa = as.integer(iwa), 
                   itask = as.integer(itask),
                   iprint = as.integer(iprint),
                   icsave = as.integer(icsave), lsave=as.logical(lsave), 
                   isave=as.integer(isave), dsave=as.double(dsave))
      itask <- result$itask
      icsave <- result$icsave
      prm <- result$x
##      cat("in lbfgsb3 parameter results:")
##      print(prm)
      g <- result$g
      iwa <- result$iwa
      wa <- result$wa
      nbd <- result$nbd
      lsave <- result$lsave
      isave <- result$isave
      dsave <- result$dsave
      if (ctrl$trace > 2) {
      cat("returned from lbfgsb3\n")
      cat("returned itask is ",itask,"\n")
      task <- tasklist[itask]
      cat("changed task to ", task,"\n")
##      task<-readline("continue")
      }

      if  (itask %in% c(4L, 20L, 21L) ) {
         if (ctrl$trace >= 2) {
          cat("computing f and g at prm=")
          print(prm)
         }
##        Compute function value f for the sample problem.
         f <- fn(prm, ...)
         nfg <- nfg + 1 # increment counter
##        Compute gradient g for the sample problem.
         if (is.null(gr)) {
             g <- numDeriv::grad(fn, prm, ...)
         } else {
             g <- gr(prm, ...)
         }
         if (ctrl$trace > 0) {
            cat("At iteration ", isave[34]," f =",f)
            if (ctrl$trace > 1) {
               cat("  max(abs(g))=",max(abs(g)))
            }
            cat("\n")
         }
         if (nfg >= ctrl$maxit) {
            if (ctrl$trace > 0) {
               cat("Exceeded function/gradient evaluation limit\n")
               break
            }
         }
      } else {
        if (itask == 1L )  { # NEW_X
##          tmp <- readline("Continue") # eventually remove this
 		##     If task is neither FG nor NEW_X we terminate execution.
          } else break
      }
 } # end repeat
## Here build return structure
##  print(result) ## only print for debugging
  info <- list(task = task, itask = itask, lsave = lsave, 
                icsave = icsave, dsave = dsave, isave = isave)
  ans <- list(prm = prm, f = f, g = g, info = info)
  ans # to return the answer visibly
##======================= The end of driver1 ============================
} # end of lbfgsb3()