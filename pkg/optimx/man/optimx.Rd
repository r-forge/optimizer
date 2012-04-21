% File src/library/stats/man/optimx.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optimx}
\alias{optimx}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optimx} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimx() does not call the SANN method.
}
\usage{
optimx(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, bdmsk=NULL, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{A vector of initial values for the parameters for which optimal 
   values are to be found. }
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.
   This includes the following methods:
		\code{"BFGS"}
		\code{"CG"}
		\code{"L-BFGS-B"}


   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }
 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{bdmsk}{Vector of indicators for bounds and masks. Masked are specified on input by
       entries in \code{bdmsk} that are 0, while free parameters have 1. On termination there
       may be values -3 for active lower bounds, and -1 for active upper bounds. 
       However, at 2011-6-28 the \code{bdmsk} vector is NOT returned. (It is returned
       by packages \code{Rvmmin} }
 \item{method}{A list of the methods to be used. 
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}.}
 \item{itnmax}{If provided as a vector of the same length as the list of methods \code{method}, 
	gives the maximum number of iterations or function values for the corresponding 
	method. If a single number is provided, this will be used for all methods. Note that
	there may be control list elements with similar functions, but this should be the
	preferred approach when using \code{optimx}.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim(). Note that 
       there are extended possibilities in \code{ucminf}, where \code{hessian} is 
       treated as an integer, but these, as well as the facility in \code{optim} are
       ignored, and the numerical hessian from the \code{kktc} package are used.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{gr}.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this. DO NOT
  use both methods. 

  Possible method codes at the time of writing are 'Nelder-Mead', 'BFGS',
  'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'Rcgmin', 'Rvmmin',
  'bobyqa', 'uobyqa', 'newuoa', 'nmkb', and 'hjkb'.

  (For convenience, 'NM' may be used instead of 'Nelder-Mead' and that is the
   display name to minimize display width. Similarly, 'LBFGSB' may be used
   instead of 'L-BFGS-B'.) 

  The default methods for unconstrained problems (no \code{lower} or
  \code{upper} specified) are 'Nelder-Mead', an implementation of the 
  Nelder and Mead (1965) and 'BFGS', a Variable Metric method based 
  on the ideas of Fletcher (1970) as modified by him in conversation 
  with J.C. Nash (see Nash, 1979). Nelder-Mead uses only function 
  values and is robust but relatively slow.  It will work reasonably 
  well for non-differentiable functions. The Variable Metric method, 
  \code{"BFGS"} updates an approximation to the inverse Hessian using 
  the BFGS update formulas, along with an acceptable point line 
  search strategy. This method appears to work best with analytic
  gradients. (\code{"Rvmmmin"} provides a box-constrained version of this
  algorithm, but there are some differences in detail so that the trajectories
  of the minimization of an unconstrained problem will not match exactly.

  If no \code{method} is given, and there are bounds constraints provided,
  the method is set to \code{"L-BFGS-B"}.

  Method \code{"CG"} is a conjugate gradients method based on that by
  Fletcher and Reeves (1964) (but with the option of Polak--Ribiere or
  Beale--Sorenson updates). The particular implementation is now dated,
  and improved yet simpler codes exist, and one of them is provided via
  \code{Rcgmin}. Furthermore this permits box constraints on parameters.
  Conjugate gradient methods will generally be more fragile than the 
  BFGS method, but as they do not store a matrix they may be successful 
  in much larger optimization problems.

  Method \code{"L-BFGS-B"} is that of Byrd \emph{et. al.} (1995) which
  allows \emph{box constraints}, that is each variable can be given a lower
  and/or upper bound. The initial value must satisfy the constraints.
  This uses a limited-memory modification of the BFGS quasi-Newton
  method. If non-trivial bounds are supplied, this method will be
  selected, with a warning.

  Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.

  The user objective function (\code{fn}) can commonly return \code{NA} 
  or \code{Inf} if the function cannot be evaluated at the supplied 
  parameters. \code{optimx} requires that the initial parameters must
  give a computable finite value of \code{fn}. However, some methods, of
  which \code{"L-BFGS-B"} is known to be a case, require that the values
  returned should always be finite. To control for excursions outside the
  domain where function \code{fn} and, if provided, \code{gr} and \code{hess},
  can be computed, \code{optimx} wraps these function using tools in 
  the package \code{optfntools}, namely \code{ufn}, \code{ugr} and \code{uhess}.

  While \code{optim} can be used recursively, and for a single parameter
  as well as many, this may not be true for \code{optimx}. \code{optim}
  also accepts a zero-length \code{par}, and just evaluates the function 
  with that argument. Such usage is strongly discouraged for \code{optimx}.

  Method \code{"nlm"} is from the package of the same name that implements
  ideas of Dennis and Schnabel (1983) and Schnabel et al. (1985). See nlm()
  for more details.

  Method \code{"nlminb"} is the package of the same name that uses the
  minimization tools of the PORT library.  The PORT documentation is at 
  <URL: http://netlib.bell-labs.com/cm/cs/cstr/153.pdf>. See nlminb()
  for details. (Though there is very little information about the methods.)

  Method \code{"spg"} is from package BB implementing a spectral projected 
  gradient method for large-scale optimization with simple constraints due
  R adaptation, with significant modifications, by Ravi Varadhan,
  Johns Hopkins University (Varadhan and Gilbert, 2009), from the original
  FORTRAN code of Birgin, Martinez, and Raydan (2001). 

  Method \code{"Rcgmin"} is from the package of that name. It implements a
  conjugate gradient algorithm with the Yuan/Dai update (ref??) and also 
  allows bounds constraints on the parameters. (Rcgmin also allows mask 
  constraints -- fixing individual parameters -- but there is no interface
  from \code{"optimx"}.) 

  Methods \code{"bobyqa"}, \code{"uobyqa"} and \code{"newuoa"} are from the 
  package \code{"minqa"} which implement optimization by quadratic approximation
  routines of the similar names due to M J D Powell (2009). See package minqa 
  for details. Note that \code{"uobyqa"} and \code{"newuoa"} are for 
  unconstrained minimization, while \code{"bobyqa"} is for box constrained
  problems.

  Methods \code{nmkb} and \code{hjkb} are box-constrained versions of the 
  Nelder-Mead implementation of Tim Kelley and the Hooke and Jeeves method. 
  When unconstrained problems are presented, functions \code{nmk} and \code{hjk}
  are called, but these naming variations are NOT permitted in the vector of
  methods \code{method} in the call to \code{optimx}. 

  The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{trace}}{Non-negative integer. If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing. trace = 0 gives no output 
      (To understand exactly what these do see the source code: higher 
      levels give more detail.)}
    \item{\code{follow.on }}{ = TRUE or FALSE. If TRUE, and there are multiple 
      methods, then the last set of 
      parameters from one method is used as the starting set for the next. }
    \item{\code{save.failures}}{ = TRUE if we wish to keep "answers" from runs 
      where the method does not return conv==0. FALSE otherwise (default).}
    \item{\code{maximize}}{ = TRUE if we want to maximize rather than minimize 
      a function. (Default FALSE). Methods nlm, nlminb, ucminf cannot maximize a
      function, so the user must explicitly minimize and carry out the adjustment
      externally. However, there is a check to avoid
      usage of these codes when maximize is TRUE. See \code{fnscale} below for 
      the method used in \code{optim} that we deprecate.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkt}}{=FALSE if we do NOT want to test the Kuhn, Karush, Tucker
      optimality conditions. The default is TRUE. However, because the Hessian
      computation may be very slow, we set \code{kkt} to be FALSE if there are 
      more than than 50 parameters when the gradient function \code{gr} is not 
      provided, and more than 500
      parameters when such a function is specified. We return logical values \code{KKT1}
      and \code{KKT2} TRUE if first and second order conditions are satisfied approximately.
      Note, however, that the tests are sensitive to scaling, and users may need
      to perform additional verification. If \code{kkt} is FALSE but \code{hessian}
      is TRUE, then \code{KKT1} is generated, but \code{KKT2} is not.}
    \item{\code{sort.result}}{= TRUE if we want to have the results table sorted 
      into reverse order of final objective function.}
    \item{\code{kkttol}}{= value to use to check for small gradient and negative
      Hessian eigenvalues. Default = .Machine$double.eps^(1/3) }
    \item{\code{kkt2tol}}{= Tolerance for eigenvalue ratio in KKT test of positive 
      definite Hessian. Default same as for kkttol }
    \item{\code{starttests}}{= TRUE if we want to run tests of the function and 
      parameters: feasibility relative to bounds, analytic vs numerical gradient, 
      scaling tests, before we try optimization methods. Default is TRUE.}
    \item{\code{dowarn}}{= TRUE if we want warnings generated by optimx. Default is 
      TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.\cr
      Default is (0.5)*.Machine$double.xmax }
    \item{\code{axsearch.tries}}{= The number of times to try the axial search. If >1, then restarts
        will be attempted until axial search does not yield a lower function value. Default is 3. 
        CAUTION: Be very careful using axial search when maximizing a function. It is much safer
        to define your own function that is (-1) * (function to be maximized) and to use minimization.}
  }

  The following \code{control} elements apply only to some of the methods. The list
  may be incomplete. See individual packages for details. 

  \describe{

    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} and \code{gr} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}. For methods from the set in
      \code{optim()}. Note potential conflicts with the control \code{maximize}.
      This also conflicts with the argument \code{ps} to the tools from 
      package optfntools.}
    \item{\code{parscale}}{A vector of scaling values for the parameters.
	Optimization is performed on \code{par/parscale} and these should be
	comparable in the sense that a unit change in any element produces
	about a unit change in the scaled value.For \code{optim}.
        This also conflicts with the argument \code{ps} to the tools from 
        package optfntools.}
    \item{\code{ndeps}}{A vector of step sizes for the finite-difference
      approximation to the gradient, on \code{par/parscale}
      scale. Defaults to \code{1e-3}. For \code{optim}.}
    \item{\code{maxit}}{The maximum number of iterations. Defaults to
      \code{100} for the derivative-based methods, and
      \code{500} for \code{"Nelder-Mead"}.}
    \item{\code{abstol}}{The absolute convergence tolerance. Only
      useful for non-negative functions, as a tolerance for reaching zero.}
    \item{\code{reltol}}{Relative convergence tolerance.  The algorithm
      stops if it is unable to reduce the value by a factor of
      \code{reltol * (abs(val) + reltol)} at a step.  Defaults to
      \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}. For \code{optim}.}
    \item{\code{alpha}, \code{beta}, \code{gamma}}{Scaling parameters
      for the \code{"Nelder-Mead"} method. \code{alpha} is the reflection
      factor (default 1.0), \code{beta} the contraction factor (0.5) and
      \code{gamma} the expansion factor (2.0).}
    \item{\code{REPORT}}{The frequency of reports for the \code{"BFGS"} and
      \code{"L-BFGS-B"} methods if \code{control$trace}
      is positive. Defaults to every 10 iterations for \code{"BFGS"} and
      \code{"L-BFGS-B"}.}
    \item{\code{type}}{for the conjugate-gradients method. Takes value
      \code{1} for the Fletcher--Reeves update, \code{2} for
      Polak--Ribiere and \code{3} for Beale--Sorenson.}
    \item{\code{lmm}}{is an integer giving the number of BFGS updates
      retained in the \code{"L-BFGS-B"} method, It defaults to \code{5}.}
    \item{\code{factr}}{controls the convergence of the \code{"L-BFGS-B"}
      method. Convergence occurs when the reduction in the objective is
      within this factor of the machine tolerance. Default is \code{1e7},
      that is a tolerance of about \code{1e-8}.}
    \item{\code{pgtol}}{helps control the convergence of the \code{"L-BFGS-B"}
      method. It is a tolerance on the projected gradient in the current
      search direction. This defaults to zero, when the check is
      suppressed.}
  }

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Note that no other attributes of \code{par}
  are copied over. (We have not verified this as at 2011-11-02.)

}
\value{

   \code{ansout <- data.frame(cbind(method=meths ,par=pars, fvalues=vals,}\cr 
   \code{fns=fevals, grs=gevals, itns=nitns, conv=convcode, KKT1=kkt1,}\cr
   \code{KKT2=kkt2, xtimes=xtimes, meths=meths))}


  A dataframe with the following columns, each having one element (row) for each
  each of the methods that has returned an answer. In each row the columne elements
  give:
  \item{method}{The name of the method used for this answer.}
  \item{par}{The best set of parameters found.}
  \item{fvalues}{The value of \code{fn} corresponding to \code{par}.}
  \item{fns}{The number of calls to \code{fn}.}
  \item{grs}{The number of calls to \code{gr}. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.}
  \item{itns}{For those methods where it is reported, the number of ``iterations''. See
    the documentation or code for particular methods for the meaning of such counts.}
  \item{conv}{An integer code. \code{0} indicates successful
    convergence. Various methods may or may not return sufficient information
	to allow all the codes to be specified. An incomplete list of codes includes
    \describe{
      \item{\code{-1}}{indicates that the function could NOT be computed at the initial
           point.} 
      \item{\code{-2}}{indicates that the bounds are inadmissible, i.e., at least one lower bound 
           exceeds the corresponding upper bound.}
      \item{\code{-3}}{indicates that at least one of the parameters is out of bounds (and the
           call has been made with control \code{allowmv2bound} = FALSE).}
      \item{\code{-4}}{indicates that the initial function returned was NOT a scalar.} 
      \item{\code{-5}}{indicates that the gradient test has failed when there is a function specified
           to compute the gradient.}
      \item{\code{-6}}{indicates that the Hessian test has failed when there is a function specified
           to compute the Hessian.}
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{3}}{indicates that the axial search failed.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
      \item{\code{10}}{indicates degeneracy of the Nelder--Mead simplex.}
      \item{\code{51}}{indicates a warning from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{52}}{indicates an error from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
    }
  }
  \item{KKT1}{A logical value returned TRUE if the solution reported has a ``small'' gradient.}
  \item{KKT2}{A logical value returned TRUE if the solution reported appears to have a 
  positive-definite Hessian.}
  \item{xtimes}{The reported execution time of the calculations for the particular method.}
  \item{meths}{A copy of the first column (method) to aid in reading a printout.}

More detail is provided in the attribute \code{"details"} to the returned answer object. 
If the returned object from optimx() is \code{ans}, this is accessed 
via the construct \cr
    \code{attr(ans, "details")} \cr
This object is a list of lists, one list for each method that has been successful 
or that has been forced by save.failures=TRUE. To get the details list for the 
third method, we use \cr
    \code{attr(ans, "details")[[3]]]]}

Each such list has possible elements:
  \item{par}{ - the final parameters of the function}
  \item{value}{ - the final value of the objective function}
  \item{convergence}{ - a code indicating whether and/or how the method terminated}
  \item{message}{ - a descriptive message, if one is provided}
  \item{fevals}{ - the number of times the objective function was evaluated during the run;
   this item has an attribute "names" that can take on the value "function"}
  \item{gevals}{ - the number of times the gradient of the objective function was evaluated 
   during the run; this item has an attribute "names" that can take on the value "gradient"}
  \item{kkt1}{ - an indicator if the first order Kuhn, Karush, Tucker condition for a 
  local minimum is satisfied. Note that we must use a tolerance to make this test, and the
  indicator is NOT definitive.}
  \item{kkt2}{  - an indicator if the second order Kuhn, Karush, Tucker condition for a 
  local minimum is satisfied. Note that we must use a tolerance to make this test, and the
  indicator is NOT definitive.}
  \item{ngatend}{ - the gradient vector computed or estimated at termination}
  \item{nhatend}{ - the Hessian matrix computed or estimated at termination}
  \item{evnhatend}{ - a vector of the eigenvalues of the Hessian at termination}
  \item{systime}{ - the system time required for execution of the method in question (named
   in "method" below). This element has the attribute "names".}
  \item{method}{ - The name of the method used. This has the attribute "CPU times (s)" which 
  itself has an attribute "names" that corresponds to the "method".}

}
\note{
  \code{optimx} will work with one-dimensional \code{par}s, but the
  default method does not work well (and will warn).  Use
  \code{\link{optimize}} instead.

  There are a series of demos available. These were set up as tests, but take quite
  a long time to execute. Once the package is loaded (via \code{require(optimx)} or
  \code{library(optimx)}, you may see available demos via 

  demo(package="optimx")

  The demo 'brown_test' may be run with the command
  demo(brown_test, package="optimx")
}
\source{
  The code for methods \code{"Nelder-Mead"}, \code{"BFGS"} and
  \code{"CG"} was based originally on Pascal code in Nash (1990) that was
  translated by \code{p2c} and then hand-optimized.  As author, I (J C Nash)
  have agreed that the code can be made freely available. However, I 
  encourage the use of some of the later methods e.g., \code{"Rcgmin"}, 
  \code{"Rvmmin"}, and \code{"nmkb"}.

  The code for method \code{"L-BFGS-B"} is based on Fortran code by Zhu,
  Byrd, Lu-Chen and Nocedal obtained from Netlib (file
  \file{opt/lbfgs\_bcm.shar}: another version is in \file{toms/778}).

  See documentation for \code{"nlm"} and \code{"nlminb"} for those
  methods, package \code{"ucminf"} for the function of the same
  name, and package BB for method \code{"spg"}.
}
\references{
  Belisle, C. J. P. (1992) Convergence theorems for a class of simulated
  annealing algorithms on \eqn{R^d}{Rd}. \emph{J Applied Probability},
  \bold{29}, 885--895.

  Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for
     convex-constrained optimization, \emph{ACM Transactions on Mathematical
     Software}. ??incomplete reference??

  Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C.  (1995) A limited
  memory algorithm for bound constrained optimization.
  \emph{SIAM J. Scientific Computing}, \bold{16}, 1190--1208.

  Dennis, J. E. and Schnabel, R. B. (1983) _Numerical Methods for
     Unconstrained Optimization and Nonlinear Equations._
     Prentice-Hall, Englewood Cliffs, NJ.

  Fletcher, R. and Reeves, C. M. (1964) Function minimization by
  conjugate gradients. \emph{Computer Journal} \bold{7}, 148--154.

  Fletcher, R (1970) A new approach to variable metric algorithms,
  \emph{Computer Journal}, \bold{13}, 317-322.

  Nash, J. C. (1979, 1990) \emph{Compact Numerical Methods for
    Computers. Linear Algebra and Function Minimisation.} Adam Hilger.

  Nelder, J. A. and Mead, R. (1965) A simplex algorithm for function
  minimization. \emph{Computer Journal} \bold{7}, 308--313.

  Nocedal, J. and Wright, S. J. (1999) \emph{Numerical Optimization}.
  Springer.

  Schnabel, R. B., Koontz, J. E. and Weiss, B. E. (1985) A modular
     system of algorithms for unconstrained minimization. _ACM Trans.
     Math. Software_, *11*, 419-440.

 Varadhan, R. and Gilbert, P.D. (2009) \bold{BB}: An R Package for Solving a Large System of
      Nonlinear Equations and for Optimizing a High-Dimensional
      Nonlinear Objective Function. \emph{Journal of Statistical Software}, 
      \bold{32}, 1--26.

}

\seealso{
  \code{\link[BB]{spg}}, \code{\link{nlm}}, \code{\link{nlminb}}, \code{\link[minqa]{bobyqa}}, 
  \code{\link[Rcgmin]{Rcgmin}}, \code{\link[Rvmmin]{Rvmmin}}, \code{\link[ucminf]{ucminf}},
  \code{\link[dfoptim]{dfoptim}}.

  \code{\link{optimize}} for one-dimensional minimization;
  \code{\link{constrOptim}} or \code{\link[BB]{spg}} for linearly constrained optimization.
}

\examples{
require(graphics)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}
cat("Rosenbrock function, with and without analytic gradients\n")

ans1<-optimx(c(-1.2,1), fr, control=list(all.methods=TRUE))
cat("Default numerical gradients via grfwd\n")
print(ans1)
# print(attr(ans1,"details"))
cat("\n\n")
cat("Analytic gradients\n")
ans2<-optimx(c(-1.2,1), fr, grr, control=list(all.methods=TRUE))
print(ans2)


cat("The next line should fail if executed because 'hessian = TRUE' no longer allowed\n")
try3<-try(ans3<-optimx(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE))
print(class(try3))

cat("\n\n")

cat("Single method call -- note control for CG is 'type=...'\n")
ans5a<-optimx(c(-1.2,1), fr, grr, method = "CG", control=list(type=1, axsearch.tries=10, trace=0))
print(ans5a)
ans5b<-optimx(c(-1.2,1), fr, grr, method = "CG", control=list(type=2, axsearch.tries=10, trace=0))
print(ans5b)
ans5c<-optimx(c(-1.2,1), fr, grr, method = "CG", control=list(type=3, axsearch.tries=10, trace=0))
print(ans5c)
cat("\n\n")


cat("25-dimensional box constrained function\n")
flb <- function(x)
    { p <- length(x); sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) }

ans7<-optimx(rep(3, 25), flb, lower=rep(2, 25), upper=rep(4, 25),
        control=list(all.methods=TRUE)) # par[24] is *not* at boundary
print(ans7)


cat("Check on maximization\n")

maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

x0<-rep(pi,4)
ans.mx<-optimx(x0,maxfn,control=list(maximize=TRUE,all.methods=TRUE,save.failures=TRUE,trace=0))
# ans.mx1<-optansout(ans.mx, filename="./ansmx.txt")
print(ans.mx)
x0b<-get.best(ans.mx)
cat("trying get.best for maxfn\n")
print(x0b)
cat("axial search on x0b\n")
## Need to be very careful with axial search when maximizing
## WRONG ## as.x0b<-axsearch(x0b, maxfn, x0b$value)
# Here is how to do it
myopcon<-optstart(length(x0b))
myopcon$MAXIMIZE=TRUE
mymax<-list(fn=maxfn, gr=NULL, hess=NULL, OPCON=myopcon)
as.x0b<-axsearch(x0b$par, ufn, -x0b$value, trace=0, fnuser=mymax)
print(as.x0b)



cat("verify using explicit negation\n")
ans.mxn<-optimx(x0,negmaxfn,control=list(all.methods=TRUE,save.failures=TRUE,trace=0))
x0bn<-get.best(ans.mxn)
print(x0bn)
as.x0bn<-axsearch(x0bn$par, negmaxfn, x0bn$value)
print(as.x0bn)
cat("\n\n")

x00<-c(1,2,3,4)
cat("Test if things work when we provide the solution!\n")
ans.mx0<-optimx(x0,maxfn,control=list(all.methods=TRUE,maximize=TRUE,save.failures=TRUE,trace=0))
optansout(ans.mx0, filename="./ansmx0.txt")
x00b<-get.best(ans.mx0)
cat("trying get.best for maxfn from solution\n")
print(x00b)
cat("\n\n")



# genrosa function code -- matches rosenbrock when npar=2 and gs=100
genrosa.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
        # Note do not at 1.0 so min at 0
	fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

genrosa.g <- function(x, gs=NULL){
# vectorized gradient for genrosa.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn1]
        # f = gs*z1*z1 + z2*z2
	gg[tn] <- 2 * (gs * z1)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 - 2 *z2 
	return(gg)
}

genrosa.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
#		z2<-1.0 - x[i-1]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}
startx<-4*seq(1:10)/3.
ans8<-optimx(startx,fn=genrosa.f,gr=genrosa.g, hess=genrosa.h, 
     control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
print(ans8)

get.result(ans8, attribute="grs")
get.result(ans8, method="spg")
cat("\n\n")


cat("Try follow-on\n")
startx<-4*seq(1:10)/3.
cat("Polyalgorithm with 200 steps NM followed by up to 75 of ucminf\n")
ans9<-optimx(startx,fn=genrosa.f,gr=genrosa.g, hess=genrosa.h, 
             method=c("Nelder-Mead","ucminf"), itnmax=c(200,75), 
             control=list(follow.on=TRUE, save.failures=TRUE,trace=0), gs=10)
print(ans9)
cat("\n\n")

cat("A case where methods do not work and we do not save failures\n")

startx<-4*seq(1:10)/3.
cat("200 steps NM is not enough to terminate\n")
ans10<-optimx(startx,fn=genrosa.f,gr=genrosa.g, method=c("Nelder-Mead"),
             itnmax=c(200), control=list(trace=0, save.failures=FALSE), gs=10)
cat("The answer should be NULL\n")
print(str(ans10))
cat("\n\n")


startx<-4*seq(1:10)/3.
cat("Try getting hessian but not kkt tests\n")
ans11<-optimx(startx,fn=genrosa.f,gr=genrosa.g, hessian=TRUE,  
    control=list(all.methods=TRUE, trace=0, save.failures=FALSE, kkt=FALSE),
    gs=10)
print(ans11)
print(attr(ans11,"details")[7])
cat("\n\n")


startx<-4*seq(1:10)/3.
cat("Use analytic hessian and no KKT tests\n")
ans12<-optimx(startx,fn=genrosa.f,gr=genrosa.g, hess=genrosa.h, 
      control=list(trace=0, save.failures=FALSE, kkt=FALSE), gs=10)
print(ans12)
print(attr(ans12,"details"))

cat("Example with bounds and masks\n")
startx<-4*seq(1:10)/7.
lo<-rep(0.5,10)
up<-rep(8,10)
bm<-c(1,1,1,1,1,0,0,1,1,1)
ans13<-optimx(startx,fn=genrosa.f,gr=genrosa.g, lower=lo, upper=up, bdmsk=bm, 
       control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=4)
print(ans13)
cat("\n\n")



# ======= check names =========
cyq.f <- function (x) {
  rv<-cyq.res(x)
  f<-sum(rv*rv)
}

cyq.res <- function (x) {
# Fletcher's chebyquad function m = n -- residuals 
   n<-length(x)
   res<-rep(0,n) # initialize
   for (i in 1:n) { #loop over resids
     rr<-0.0
     for (k in 1:n) {
	z7<-1.0
	z2<-2.0*x[k]-1.0
        z8<-z2
        j<-1
        while (j<i) {
            z6<-z7
            z7<-z8
            z8<-2*z2*z7-z6 # recurrence to compute Chebyshev polynomial
            j<-j+1
        } # end recurrence loop
        rr<-rr+z8
      } # end loop on k
      rr<-rr/n
      if (2*trunc(i/2) == i) { rr <- rr + 1.0/(i*i - 1) }
      res[i]<-rr
    } # end loop on i
    res
}

cyq.jac<- function (x) {
#  Chebyquad Jacobian matrix
   n<-length(x)
   cj<-matrix(0.0, n, n)
   for (i in 1:n) { # loop over rows
     for (k in 1:n) { # loop over columns (parameters)
       z5<-0.0
       cj[i,k]<-2.0
       z8<-2.0*x[k]-1.0 
       z2<-z8
       z7<-1.0
       j<- 1
       while (j<i) { # recurrence loop
         z4<-z5
         z5<-cj[i,k]
         cj[i,k]<-4.0*z8+2.0*z2*z5-z4
         z6<-z7
         z7<-z8
         z8<-2.0*z2*z7-z6
         j<- j+1
       } # end recurrence loop
       cj[i,k]<-cj[i,k]/n
     } # end loop on k
   } # end loop on i
   cj
}


cyq.g <- function (x) {
   cj<-cyq.jac(x)
   rv<-cyq.res(x)
   # Need to escape the \% in .Rd file
   gg<- as.vector(2.0* rv \%*\% cj)
}

startx<-c(0.2, 0.6)
names(startx)<-c("First","Second")
cat("check the names on the parameters\n")
print(startx)
ansnames<-optimx(startx, cyq.f, cyq.g, control=list(all.methods=TRUE))
ansnames
cat("check the names on the parameters\n")
print(ansnames$method)
print(ansnames$par)
bestans<-get.best(ansnames)
cat("The get.best result\n")
print(bestans)
cat("\n\n")

startx<-c(0.2, 0.4, 0.6, 0.8)
cat("check Hessian -- no bounds\n")
# should check bounds, control.kkt=FALSE, etc.
anshes<-optimx(startx,cyq.f, cyq.g, hessian=TRUE, control=list(all.methods=TRUE))
# hesBFGS<-get.result(anshes, attribute="nhatend")
adet<-attr(anshes,"details")
nmeth<-length(adet)
for (i in 1:nmeth){
   print(adet[[i]]$method)
   print(adet[[i]]$nhatend)
   print(adet[[i]]$evnhatend)
   cat("=======================================\n")
}

}
\keyword{nonlinear}
\keyword{optimize}
