% File src/library/stats/man/optim.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optimr}
\alias{optimr}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optim} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimz::optim() does not call the SANN method.
}
\usage{
optimr(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{method}{A list of the methods to be used. 
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}.
       The default of NULL causes an appropriate set of methods to be supplied
       depending on the presence or absence of bounds on the parameters. The default
       unconstrained set is Rvmminu, Rcgminu, lbfgsb3, newuoa and nmkb.
       The default bounds constrained set is Rvmminb, Rcgminb, lbfgsb3, bobyqa and nmkb.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim().}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this. DO NOT
  use both methods. 

       "Rtnmin", "hjn", "lbfgs")

  Possible method codes are 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 
  'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 
  'lbfgsb3', 'Rtnmin', 'hjn', 'lbfgs' or 'Rvmmin'. These are in base R or in CRAN repositories. 
  From R-forge, method 'Rtnmin' is available. Other methods are likely to
  be added over time.

  The default methods for unconstrained problems (no \code{lower} or
  \code{upper} specified) are an implementation of the Nelder and Mead
  (1965) and a Variable Metric method based on the ideas of Fletcher
  (1970) as modified by him in conversation with Nash (1979). Nelder-Mead
  uses only function values and is robust but relatively slow.  It will 
  work reasonably well for non-differentiable functions. The Variable
  Metric method, \code{"BFGS"} updates an approximation to the inverse
  Hessian using the BFGS update formulas, along with an acceptable point
  line search strategy. This method appears to work best with analytic
  gradients. (\code{"Rvmmmin"} provides a box-constrained version of this
  algorithm.

  If no \code{method} is given, and there are bounds constraints provided,
  the method is set to \code{"L-BFGS-B"}.

  Method \code{"CG"} is a conjugate gradients method based on that by
  Fletcher and Reeves (1964) (but with the option of Polak--Ribiere or
  Beale--Sorenson updates). The particular implementation is now dated,
  and improved yet simpler codes have been implemented. Furthermore, 
  \code{"Rcgmin"} allows box constraints as well as fixed (masked)
  parameters. Conjugate gradient methods will generally be more fragile 
  than the BFGS method, but as they do not store a matrix they may be 
  successful in optimization problems with a large number of parameters.

  Method \code{"L-BFGS-B"} is that of Byrd \emph{et. al.} (1995) which
  allows \emph{box constraints}, that is each variable can be given a lower
  and/or upper bound. The initial value must satisfy the constraints.
  This uses a limited-memory modification of the BFGS quasi-Newton
  method. If non-trivial bounds are supplied, this method is selected
  by the original \code{optim()} function, with a warning. Unfortunately,
  the authors of the original Fortran version of this method released a
  correction for bugs in 2011, but these have not been incorporated into
  the distributed R codes, which are a C translation of a version that
  appears to be from the mid-1990s. While it seems the errors affect 
  very few computations, users may wish to use the Fortran codes in 
  package \code{lbfgsb3}. 

  Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.

  Function \code{fn} can return \code{NA} or \code{Inf} if the function
  cannot be evaluated at the supplied value, but the initial value must
  have a computable finite value of \code{fn}. However, some methods, of
  which \code{"L-BFGS-B"} is known to be a case, require that the values
  returned should always be finite.

  While \code{optim} can be used recursively, and for a single parameter
  as well as many, this may not be true for \code{optimx}. \code{optim}
  also accepts a zero-length \code{par}, and just evaluates the function 
  with that argument.

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
  constraints -- fixing individual parameters -- but there is as yet no 
  interface from \code{"optimr"}.) 

  Method \code{"Rvmmin"} is from the package of that name. It implements 
  the same variable metric method as the base optim() function with method
  \code{"BFGS"} but allows bounds constraints on the parameters. (Rvmmin 
  also allows mask constraints -- fixing individual parameters -- but 
  there is as yet no interface from \code{"optimr"}.) 

  Method \code{"Rtnmin"} is from the package of that name. It implements a
  truncated Newton method of Stephen Nash translated from Matlab. It 
  allows bounds constraints on the parameters. 

  Methods \code{"bobyqa"}, \code{"uobyqa"} and \code{"newuoa"} are from the 
  package \code{"minqa"} which implement optimization by quadratic approximation
  routines of the similar names due to M J D Powell (2009). See package minqa 
  for details. Note that \code{"uobyqa"} and \code{"newuoa"} are for 
  unconstrained minimization, while \code{"bobyqa"} is for box constrained
  problems. While \code{"uobyqa"} may be specified, it is NOT part of the 
  \code{all.methods = TRUE} set.

  Methods \code{"nmkb"} and \code{"hjkb"} are from package dfoptim. They 
  implement respectively variants of the Nelder-Mead and Hooke and Jeeves
  derivative-free methods, but both allow bounds constraints. However, it is
  important to note that \code{"nmkb"} must NOT have starting parameters on
  a lower or upper bound, as a transformation of the paramters is used to
  effect the constraints. 

  Method \code{ucminf} is an unconstrained minimizer based on a
     BFGS update of the inverse Hessian approximation, much like
     that of methods 'BFGS' or 'Rvmmin' but with a more 
     sophisticated line search. See Nielsen (2000).

  Method \code{lbfgs} is an unconstrained minimizer based on 
     the libLBFGS library by Naoaki Okazaki, in turn based on an
     implementation of the L-BFGS method written by Jorge Nocedal.


  Method \code{"hjn"} is a conservative implementation of a Hooke and 
  Jeeves (1961) 

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
      where the method does not return convcode==0. FALSE otherwise (default).}
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
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkttol}}{= value to use to check for small gradient and negative
      Hessian eigenvalues. Default = .Machine$double.eps^(1/3) }
    \item{\code{kkt2tol}}{= Tolerance for eigenvalue ratio in KKT test of positive 
      definite Hessian. Default same as for kkttol }
    \item{\code{starttests}}{= TRUE if we want to run tests of the function and 
      parameters: feasibility relative to bounds, analytic vs numerical gradient, 
      scaling tests, before we try optimization methods. Default is TRUE.}
    \item{\code{dowarn}}{= TRUE if we want warnings generated by optimx. Default is 
      TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.
      Default is (0.5)*.Machine$double.xmax }
    \item{\code{usenumDeriv}}{= TRUE if the \code{numDeriv} function \code{grad()} is
      to be used to compute gradients when the argument \code{gr} is NULL or not supplied.}
  }

  The following \code{control} elements apply only to some of the methods. The list
  may be incomplete. See individual packages for details. 

  \describe{

    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} and \code{gr} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}. For methods from the set in
      \code{optim()}. Note potential conflicts with the control \code{maximize}.}
    \item{\code{parscale}}{A vector of scaling values for the parameters.
	Optimization is performed on \code{par/parscale} and these should be
	comparable in the sense that a unit change in any element produces
	about a unit change in the scaled value.For \code{optim}.}
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
  are copied over. (We have not verified this as at 2009-07-29.)

}
\value{

     For ‘optim’, a list with components:

   \item{par}{The best set of parameters found.}
   \item{value}{The value of ‘fn’ corresponding to ‘par’.}
   \item{counts}{ A two-element integer vector giving the number of calls to
          ‘fn’ and ‘gr’ respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to ‘fn’
          to compute a finite-difference approximation to the gradient.}
   \item{convergence}{ An integer code. ‘0’ indicates successful completion}
   \item{ message}{ A character string giving any additional information returned
          by the optimizer, or ‘NULL’.}
   \item{hessian}{ Always NULL for this routine.}
}
\source{
See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.
}
\references{

 See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

Dai YH, and Yuan Y (2001). An efficient hybrid conjugate 
      gradient method for unconstrained optimization. Annals of 
       Operations Research 103 (1-4), 33–47.

 Hooke R. and Jeeves, TA (1961). Direct search solution of numerical and statistical problems. 
   Journal of the Association for Computing Machinery (ACM). 8 (2): 212–229.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

 Nielsen HB (2000). UCMINF - An Algorithm For Unconstrained,
     Nonlinear Optimization, Report IMM-REP-2000-18, Department of
     Mathematical Modelling, Technical University of Denmark. <URL:
     http://www2.imm.dtu.dk/~hbn/publ/TR0019.ps> or <URL:
     http://orbit.dtu.dk/recid/200975>.

 Nocedal J, and Wright SJ (1999). Numerical optimization. New York: Springer. 2nd Edition 2006.
 
}

\keyword{nonlinear}
\keyword{optimize}
