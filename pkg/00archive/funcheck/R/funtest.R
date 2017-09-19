funtest <- function(xpar, fn=NULL, gr=NULL, hess=NULL, res=NULL, jac=NULL, rsd=NULL, lower=NULL, 
                    upper=NULL, cctrl=list(trace=1), ... )
{
##funtest <- function(xpar, fn, gr, hess, res, jac, rsd, lower=rep(-Inf,length(x)), upper=rep(Inf,length(x)), 
## bdmsk=rep(1,length(xpar)), ... )
#
#  A function to test the function (fn) and gradient (gr) to be used in the optimx or related packages in R
# ?? NOT TRUE!  This function can take-in multiple starting values
#
# Input:
#  xpar = a vector of starting values
#  fn = objective function (assumed to be sufficeintly differentiable)
#  gr = function to compute gradient of objective function
#  hessian = function to compute hessian of objective function
#  res = gradient of objective function
#  jac = gradient of objective function
#  rsd = gradient of objective function
#  lower, upper = Lower and upper bounds on the variables for the problem at hand. 
#                 Not all methods can handle bounds.
#  bdmsk = an indicator vector of which variables (in par) are fixed or "masked", that is, are not allowed
#          to change during the optimisation
#
# Output:
#  ?? what will the function return??
#
#  Authors:  Ravi Varadhan & John Nash
#  Date:  June 19, 2009
#################################################################
  maxard10<-function(one, two) { # get maximum absolute relative difference scaled by 10. in denominator
# This internal function is used to make comparisons using a relative difference, but avoiding zero divide
    result<-max(abs((one-two)/(abs(one)+abs(two)+10.0)))
    return(result)
  }
########## end internal function(s) ############
  require("numDeriv")  ## make sure numerical derivatives available
  x<-xpar # ?? clean up later -- just for compatibility
  answer<-list()
# make sure we inform user when things not provided
  if (is.null(fn)) answer$fn<-NA
  if (is.null(gr)) answer$gr<-NA
  if (is.null(hess)) answer$hess<-NA
  if (is.null(res)) answer$res<-NA
  if (is.null(jac)) answer$jac<-NA
  if (is.null(rsd)) answer$rsd<-NA

  xargs<-list(...)
  answer$xargs<-xargs

  if(cctrl$trace>0) {
  # ?? Do we need sxargs?
    sxargs<-paste(names(xargs)[1],"=",unlist(xargs)[1],sep='')
    if (length(names(xargs))>1) {
       for (j in 2:length(names(xargs))) {
          sxargs<-paste(sxargs,",",names(xargs)[j],"=",unlist(xargs)[j],sep='')
       }
    }
    cat("function extra (dot-dot-dot) arguments:",sxargs,"\n")
  }
  # parameter vectors
  # Check parameters are in right form
  if(is.null(x)) stop("Null parameter vector")
  if(!is.null(dim(x))) stop("Parameter should be a vector, not a matrix!", call. = FALSE)
  if (! is.vector(x) ) {
	stop("The parameters are NOT in a vector")
  }
  answer$par<-x
  n<-length(x)
  answer$n<-n
  if (cctrl$trace>0) cat("Function has ",n," parameters\n")
  answer$lower<-NA
  answer$upper<-NA
  if (cctrl$trace > 0) {
     if (! is.null(lower) ) {
         cat("Lower bounds:")
         print(lower)
     }
     if (! is.null(upper) ) {
         cat("upper bounds:")
         print(upper)
     }
  }
  if (is.null(lower) && is.null(upper)) nobounds<-TRUE
  answer$lower<-lower
  if (is.null(lower)) {
     if (cctrl$trace > 0) cat("No lower bounds\n")
     lower<-rep(-Inf,n)
  }
  answer$upper<-upper
  if (is.null(upper)) {
     if (cctrl$trace > 0) cat("No upper bounds\n")
     upper<-rep(Inf,n)
  }
  if (cctrl$trace > 0) cat("Initial parameters:")
  if (cctrl$trace > 0) print(x)
#  if (cctrl$trace > 0) cat("bdmsk:") # Improve this ???
#  if (cctrl$trace > 0) print(bdmsk)
  ## the dim() of a scalar is NULL, so scalar starting value gives TRUE in if() test.
  ## We used this when we tried a matrix of input parameters i.e, rows are each a vector of starting values
#?? CODE HERE TO BE COMMON WITH optimx() and funtst()   090618 -- but it is not at present
#    answer$bdmsk<-bdmsk # ?? do we need to ensure integer??
  infeasible<-as.logical(FALSE)
  bstate<-vector(mode="character", length=n)
  for (i in 1:n) {
    if ( (lower[i]<=x[i]) && (x[i]<=upper[i])) {
      bstate[i]<-" In Bounds "
    } else { 
#     if(bdmsk[i]!=0) {
        infeasible<-TRUE
#     } # Note that masked parameters may be out of bounds -- should inform user
      if (lower[i]>x[i]) {bstate[i]<-" Out of Bounds LOW" } else { bstate[i]<-" Out of Bounds HIGH " }
    } # end if in bounds
#    if (cctrl$trace > 0) cat("x[",i,"]: ",lower[i],"  <?",x[i],"  <?",upper[i],"  ",bdmsk[i],"   ",bstate,"\n")
    if (cctrl$trace > 0) cat(lower[i],"  <?","x[",i,"] =  ",x[i],"  <?",upper[i],"  ",bstate[i],"\n")
  } # end of for loop over parameter vector elements
  answer$bstate<-bstate
  if(infeasible) { 
    stop("Infeasible point, no further tests") # for moment, just stop, but later need better exit mechanism
  } # ?? else ...

  # arguments to call
  # Feasible point, check function
  fval<-try(fn(x,...)) # ?? we need to sort out how to do the ... args right
  # Check if function can be computed
  if (class(fval) == "try-error") {
    infeasible <- TRUE
    print(fval)
    stop("Cannot evaluate function at initial parameters")
  }
  # Also check that it is returned as a scalar
  if ((is.vector(fval) && (length(fval)>1)) || is.list(fval) || is.matrix(fval) || is.array(fval) || ! is.numeric(fval) ) {
      stop("Function provided is not returning a scalar number")
  }
  
  if (is.infinite(fval)) {
     stop("Function returned is infinite (non-computable)")
  }
  
  if (cctrl$trace>0) cat("function value = ",fval,"\n")
  answer$fval<-fval

  answer$grd<-NA # set NA (??NULL), then change when found
  # Test if gradient code present
  if(! is.null(gr)) {
    gn<-grad(fn, x,...)
    ga<-gr(x,...)
    grd<-maxard10(gn, ga)
    if (cctrl$trace > 0) cat("Max rel diff shift 10 for gradient = ",grd,"\n")
    answer$grd<-grd
  } else {
     if(cctrl$trace >0) cat("Gradient not supplied \n")
  }
  # Test if residuals present
  answer$Jrd<-NA
  answer$rsdrd<-NA
  if(is.null(res)) {
    msg<-paste("No residual code, so no Jacobian tests",sep='')
    if(cctrl$trace>0) cat(msg,"\n")
  } else {
  # If Jacobian code is present, test code
    if(! is.null(jac)) {
      Jn<-jacobian(res,x,...)
      m<-nrow(Jn) # number of residuals / rows in Jacobian
      Ja<-jac(x,...)
      Jrd<-maxard10(Jn, Ja)
      if (cctrl$trace > 0) cat("Max rel diff shift 10 for Jacobian = ",Jrd,"\n")
      answer$Jrd<-Jrd        # If residual second derivatives present, test code    
      if(! is.null(rsd)) {
         rsda<-rsd(x,...)
         rsdn<-array(0, c(m,n,n)) # set up array for numerical residual second derivatives
         # Create rsdn locally -- may want to adjust things
         hh<-1e-7 # Choose a fixed shift for derivs. Not necessarily a good idea??
         for (j in 1:n) {
            xadj<-x
            xadj[j]<-xadj[j]+hh
            Jj<-jac(xadj,...)
            rsdn[,,j]<-(Jj-Ja)/hh
         }
         rsdrd<-maxard10(rsdn, rsda)
         if (cctrl$trace > 0) cat("Max rel diff shift 10 for residual second derivatives = ",rsdrd,"\n")
         answer$rsdrd<-rsdrd
      } else {# end residual second derivative
         if (cctrl$trace > 0) cat("No residual second derivative code provided\n")
      }
    } else { # end jacobian stuff
      if (cctrl$trace > 0) cat("No Jacobian code provided\n")
    }
  } # end of residuals
  # Test hessian if present
  answer$Hrd<-NA
  if( is.null(hess)) {
    if (cctrl$trace > 0) cat("No Hessian code provided\n")
  } else {
     Hn<-hessian(fn, x, ...)
     Ha<-hess(x,...)
     Hrd<-maxard10(Hn, Ha)
     if(cctrl$trace > 0) cat("Max rel diff shift 10 for Hessian = ",Hrd,"\n")
     answer$Hrd<-Hrd
  }
  return(answer) # ?? have not decided all content
}
### end of funtest ***

