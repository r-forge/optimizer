optest.check <- function(par, ufn, ugr=NULL, uhess=NULL, lower=-Inf, upper=Inf,  
                             ctrl=NULL, ...) {

## Should be run whenever we are not sure parameters and function are
## admissible. 
##
## Inputs: ?? 
##   control list -- ctrl. ufn, ugr

## Outputs: ?? failed-checks info.

###############################################################################


# Check parameters are in right form
  if (!is.null(dim(par))) stop("Parameters should be a vector, not a matrix!")
  if (! is.vector(par) ) {
	stop("The parameters are NOT in a vector")
  }
  npar<-length(par)
  if (is.null(ctrl)) ctrl <- ctrldefault(npar)
  optchk<-list() # set up output kust of the checks
#  if (ctrl$starttests) {
     # Check parameters in bounds (090601: As yet not dealing with masks ??)
     infeasible<-FALSE
     if (ctrl$trace > 1) cat("Function has ",npar," arguments\n")
     if (! ctrl$have.bounds) { # Don't do the check if we already know there are bounds
          
        shift2bound <- ! ctrl$keepinputpar
        bc <- bmchk(par, lower=lower, upper=upper, trace=ctrl$trace, shift2bound)
        if (! bc$admissible) stop("At least one lower bound is > corresponding upper bound")
        if (infeasible && ctrl$dowarn) warning("Parameters may be out of bounds")
        if (ctrl$trace > 0) {
           cat("Parameter relation to bounds\n")
           print(bc$bchar)
        }
        if (bc$parchanged) par <- bc$bvec
     }

     # Check if function can be computed
     checkfn <- fnchk(par, ufn, trace=ctrl$trace, ...)
     if (checkfn$infeasible) {
        cat("fnchk exit code and msg:",checkfn$excode," ",checkfn$msg,"\n")
        stop("Cannot evaluate function at initial parameters")
     }

     optchk$grbad <- FALSE
     if (! is.null(ugr) && ! ctrl$usenumDeriv && ! is.character(ugr)){ # check gradient
       gname <- deparse(substitute(ugr))
       if (ctrl$trace > 0) cat("Analytic gradient from function ",gname,"\n\n")
          fval <- ufn(par,...) 
          gn <- grad(func=ufn, x=par,...) # 
          ga <- ugr(par, ...)
#130929          badgrad<-TRUE
#130929          if (all(! is.na(ga)) & all(is.finite(ga))) badgrad<-FALSE
          # Now test for equality (090612: ?? There may be better choices for the tolerances.
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(gn-ga))/(1 + abs(fval)) >= teps) {
            stop("Gradient function might be wrong - check it! \n", call.=FALSE)
            optchk$grbad <- TRUE # Never get here if we stop ??
          }
       } else if (ctrl$trace > 0) cat("Analytic gradient not made available.\n")

       optchk$hessbad <- FALSE
##?? need to fix this
       if (! is.null(uhess) && ! is.character(uhess)){ # check Hessian - if character then numeric
          hname <- deparse(substitute(uhess))
          if (ctrl$trace > 0) cat("Analytic hessian from function ",hname,"\n\n")
          hn <- hessian(func=ufn, x=par,...) # ?? should we use dotdat
          ha <- uhess(par, ...)
          # Now test for equality
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
          optchk$hessbad <- TRUE
       } else if (ctrl$trace > 0) cat("Analytic Hessian not made available.\n")
   
# Scaling check  091219
    if (ctrl$starttests) {
        optchk$scalebad <- FALSE
	srat<-scalecheck(par, lower, upper,ctrl$dowarn)
	sratv<-c(srat$lpratio, srat$lbratio)
	if (max(sratv,na.rm=TRUE) > ctrl$scaletol) { 
		warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
		if (ctrl$dowarn) warning(warnstr)
             optchk$scalebad <- TRUE
	}
        if (ctrl$trace > 0) {
		cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
	}
    }
# end scaling check

## ?? what to return
    optchk
} ## end of optimx.check

