optimx.check <- function(par, ufn, ugr, uhess, lower=-Inf, upper=Inf, 
             hessian=FALSE, ctrl, have.bounds=FALSE, usenumDeriv=FALSE, ...) {
## Should be run whenever we are not sure parameters and function are
## admissible. 
##
## Inputs:  
# par - a vector of initial values for the parameters 
# ufn  - A function to be minimized (or maximized)
# ugr  - A function to return (as a vector) the gradient 
# uhess- A function to return (as a symmetric matrix) the Hessian of the objective 
# lower, upper - Bounds on the variables
# control - A list of control parameters. 
# ... - further arguments to be passed to fn and gr

## Output: optchk -- a list of three elements
#     - grbad -- TRUE if gradient failed-check
#     - hessbad -- TRUE if hessian failed-check
#     - scalebad -- TRUE if scale of parameters or bounds too severe

###############################################################################
## Code more or less common to funtest, funcheck and optimx <<<
# Check parameters are in right form
  if (!is.null(dim(par))) stop("Parameter should be a vector, not a matrix!", call. = FALSE)
  if (! is.vector(par) ) {
	stop("The parameters are NOT in a vector")
  }
  npar<-length(par)
  optchk<-list(grbad=FALSE, hessbad=FALSE, scalebad=FALSE) # output of the checks
  if (ctrl$starttests) {
	# Check parameters in bounds (090601: As yet not dealing with masks?)
	#  bdmsk<-as.vector(bdmset[k, ])
	infeasible<-FALSE
	if (ctrl$trace > 0) cat("Function has ",npar," arguments\n")
	if (have.bounds) {
    	  # Expand bounds to vectors if needed
          # Note 20100610: we do not check if there is a vector of wrong length.!!?
    	  if (length(lower)==1 ) lower <- rep(lower, npar)
    	  if (length(upper)==1 ) upper <- rep(upper, npar)
    	  bstate<-vector(mode="character", length=npar)
    	  for (i in 1:npar) {
     	    if ( (lower[i]<=par[i]) && (par[i]<=upper[i])) {
    	      bstate[i]<-" In Bounds "
            } else { 
            #   if (bdmsk[i]!=0) {
                  infeasible<-TRUE
            #   }
                if (lower[i]>par[i]) {bstate[i]<-" Out of Bounds LOW" } else { bstate[i]<-" Out of Bounds HIGH " }
            } # end if in bounds
            if (ctrl$trace > 0) cat("par[",i,"]: ",lower[i],"  <?",par[i],"  <?",upper[i],"  ",bstate[i],"\n") # fix to add index 150604
          } # end of for loop over parameter vector elements
	  if (infeasible) { ## maybe don't want to stop!!?
        	stop("Infeasible point, no further tests")
	  } 
  	} # end have.bounds
        # Check if function can be computed
        firsttry<-try(finit<-ufn(par), silent=TRUE ) # 20100711
        # Note: This incurs one EXTRA function evaluation because optimx is a wrapper for other methods
        if (inherits(firsttry, "try-error")) {
    	   infeasible <- TRUE
           stop("Cannot evaluate function at initial parameters")
        }
        # Also check that it is returned as a scalar
       if (!(is.vector(finit) && (length(finit)==1)) || is.list(finit) || 
            is.matrix(finit) || is.array(finit) || ! is.numeric(finit) ) {
           stop("Function provided is not returning a scalar number")
       }
       if (is.infinite(finit) || is.na(finit)) {
          stop("Function returned is infinite or NA (non-computable)")
       }
  }

  if (ctrl$starttests && ! is.null(ugr)) { # add check to see if ugr present
     if (! is.null(ugr) && ! usenumDeriv && ! is.character(ugr)){ # check gradient
       gname <- deparse(substitute(ugr))
       if (ctrl$trace>0) cat("Analytic gradient from function ",gname,"\n\n")
          fval <- ufn(par) 
          gn <- numDeriv::grad(func=ufn, x=par) # 211015 Is this the problem? CHANGED
          ga <- ugr(par)
          # Now test for equality (090612: There may be better choices for the tolerances.
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(gn-ga))/(1 + abs(fval)) >= teps) {
            # stop("Gradient function might be wrong - check it! \n", call.=FALSE)
            optchk$grbad <- TRUE # Never get here if we stop 
          }
       } else if (ctrl$trace>0) cat("Analytic gradient not made available.\n")

       if (! is.null(uhess) && ! is.character(uhess)){ # check Hessian - if character then numeric
          hname <- deparse(substitute(uhess))
          if (ctrl$trace>0) cat("Analytic hessian from function ",hname,"\n\n")
          hn <- numDeriv::hessian(func=ufn, x=par) # dotargs are in ufn
          ha <- uhess(par)
          # Now test for equality
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) {
             # stop("Hessian function might be wrong - check it! \n", call.=FALSE)
             optchk$hessbad <- TRUE
          }
       } else if (ctrl$trace>0) cat("Analytic Hessian not made available.\n")
   }
# Scaling check  091219
    if (ctrl$starttests) {
	srat<-scalecheck(par, lower, upper,ctrl$dowarn)
	sratv<-c(srat$lpratio, srat$lbratio)
	if (max(sratv,na.rm=TRUE) > ctrl$scaletol) { 
		warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
		if (ctrl$dowarn) warning(warnstr)
             optchk$scalebad <- TRUE
	}
        if (ctrl$trace>0) {
		cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
	}
    }
# end scaling check
## return
    optchk
} ## end of optimx.check

