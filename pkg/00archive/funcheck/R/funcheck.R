funcheck <- function(xpar, fname, lower=NULL, upper=NULL, cctrl=list(trace=1), ... )
{
# funcheck <- function(xpar, fname, lower=NULL, upoer=NULL, bdmsk=NULL, cctrl=list(trace=1), ... )
#  A function to check the nonlinear optimization file that is "fname.R", where fname
#  is provided by the user.
#  The intention is to automatically test the gradient, hessian, Jacobian, Jacobian second derivatives,
#    as well as bounds
#
#  This function can take-in multiple starting values
#
# Input:
#  xpar = a vector of starting values
#      If NULL, tries to use fname.setup to generate them.
#  fname = objective function (assumed to be sufficeintly differentiable)
#  lower = vector giving lower bounds on parameters in xpar or NULL
#     Note: if xpar is NULL, lower bounds are taken from the function identified by fname.setup
#       if these are set.
#  upper = vector giving upper bounds on parameters in xpar or NULL
#     Note: if xpar is NULL, upper bounds are taken from the function identified by fname.setup
#       if these are set.
#  ?? at the moment bdmsk is NOT used
#  cctrl = a list of control information FOR THE CHECKING PROGRAM. See Details.
#          The name has been changed from control to avoid confusion with control list in optim/optimx
#  ...     = other arguments to the function identified by fname
#
#
# Output:
#  ?? what will the function return??
#
#  Authors:  Ravi Varadhan & John Nash
#  Date: June 18, 2009
#################################################################
  maxard10<-function(one, two) { # get maximum absolute relative difference scaled by 10. in denominator
# This internal function is used to make comparisons using a relative difference, but avoiding zero divide
    result<-max(abs((one-two)/(abs(one)+abs(two)+10.0)))
    return(result)
  }
########## end internal function(s) ############
  require("numDeriv")  ## make sure numerical derivatives available
  pfile<-paste(fname,".R",sep='')
  if (cctrl$trace > 0) cat("pfile:",pfile,"\n")
  dtxt<-dir(pattern=pfile)
## NOTE: Only source() loaded files are checked here -- What about library ones??
  if (any(pfile == dtxt)) { # found file 
      source(pfile)
      if (cctrl$trace > 0) cat("Function ",pfile," has been (re-)loaded\n")
  } else {
     msg<-paste("Function file ",pfile," does not appear to be present here",sep='')
     stop(msg)
  }
  answer<-list()
  answer$fname<-fname
  xargs<-list(...)
  answer$xargs<-xargs
  sxargs<-paste(names(xargs)[1],"=",unlist(xargs)[1],sep='')
  if (length(names(xargs))>1) {
     for (j in 2:length(names(xargs))) {
        sxargs<-paste(sxargs,",",names(xargs)[j],"=",unlist(xargs)[j],sep='')
     }
  }

  if(cctrl$trace>0) {
      cat("function extra (dot-dot-dot) arguments:",sxargs,"\n")
  }
  # parameter vectors
  if (is.null(xpar)) { 
     if (cctrl$trace > 0) cat("Using function setup for data for tests\n")
     funs<-paste(fname,".setup",sep='')
     if (cctrl$trace > 0) cat("Using ",funs,"() for setup\n", sep="")
     lstxt <- ls(pattern=funs, name=".GlobalEnv")
     if (cctrl$trace > 0) print(lstxt)
     if (cctrl$trace > 0) cat("ls gives: ",lstxt,"\n")
     if( any(funs == lstxt) ) {
       if (cctrl$trace > 0) cat("found ",funs," loaded\n")
       setup<-try(eval(call(funs, fargs))) # ?? need to nicely incorporate the dot-dot-dot == fargs variables??
       if (cctrl$trace > 0) {
         cat("Result of ",funs,"\n")
         print(setup)
         cat("==============================\n")
       }
       x<-setup$x
       if (! is.null(lower) || ! is.null(upper)) {
           msg<-paste("Bounds supplied in function call are ignored.",sep='')
           if (cctrl$trace>0) cat(msg,"\n")
           warning(msg)
#       ?? put msg into answer??
       }
       lower<-setup$lower # Note: we ignore any bounds in funcheck call.
       upper<-setup$upper
#       bdmsk<-setup$bdmsk
       fargs<-setup$fargs  # these are the dot-dot-dot variables
     } else {
       msg<-paste("Function ",funs," not found in Global environment",sep='')
       stop(msg)
     }
  } else {
     if (cctrl$trace > 0) cat("Using data inputs for tests\n")
     x<-xpar # copy the values
  }
# Check parameters are in right form
  answer$par<-x
  if(!is.null(dim(x))) stop("Parameter should be a vector, not a matrix!", call. = FALSE)
  if (! is.vector(x) ) {
	stop("The parameters are NOT in a vector")
  }
  n<-length(x)
  answer$n<-n
  answer$lower<-NA
  answer$upper<-NA
  if (cctrl$trace>0) cat("Function has ",n," parameters\n")
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

  # Feasible point, check function


  funf<-paste(fname,".f",sep='')
  answer$fval<-NA # set NA (??NULL), then change when found
  # Test if function code present
  if(any(funf==ls(pattern=funf, name=".GlobalEnv"))) {
    fcall<-paste(funf,"(x,...)",sep='')
  #  cat("fcall:",fcall,"\n")
    if (cctrl$trace > 0) cat("about to call function ",funf,"\n")
    fval<-try(eval(parse(text=fcall)))  # Need eval because funf is character string
    # Check if function can be computed
    if (class(fval) == "try-error") {
       infeasible <- TRUE
       print(fval)
       stop("Cannot evaluate function at initial parameters")
    }
  } else {
    stop("Function code not found -- will not continue\n") # ?? eventually do so -- may be res and jac??
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
  
  fung<-paste(fname,".g",sep='')
  answer$grd<-NA # set NA (??NULL), then change when found
  # Test if gradient code present
  if(any(fung==ls(pattern=fung, name=".GlobalEnv"))) {
    if (cctrl$trace > 0) cat("about to call grad from numDeriv with ",funf,"\n")
    fcall<-paste("grad(",funf,",x,...)",sep='')
#    cat("fcall:",fcall,"\n")
    gn<-try(eval(parse(text=fcall)))  # Need eval because funf is character string
    # cat("gn:")
    # print(gn)
    fcall<-paste(fung,"(x,...)",sep='') # ?? handle ... better??
    ga<-try(eval(parse(text=fcall))) # ?? dots vars?
    grd<-maxard10(gn, ga)
    if (cctrl$trace > 0) cat("Max rel diff shift 10 for gradient = ",grd,"\n")
    answer$grd<-grd
  } else {
     msg<-paste("Function ",funf," not found in global environment",sep='')
     if(cctrl$trace >0) cat("Warning: ",msg,"\n")
  }
  # Test if residuals present
  funres<-paste(fname,".res",sep='')
  answer$Jrd<-NA
  answer$rsdrd<-NA
  if(! any(funres==ls(pattern=funres, name=".GlobalEnv"))) {
     msg<-paste("No residual code, so no Jacobian tests",sep='')
     if(cctrl$trace>0) cat(msg,"\n")
  } else {
  # If Jacobian code is present, test code
    funjac<-paste(fname,".jac",sep='')
    if(any(funjac==ls(pattern=funjac, name=".GlobalEnv"))) {
      fcall<-paste("jacobian(",funres,",x,...)",sep='')
#      cat("fcall:",fcall,"\n")
      Jn<-try(eval(parse(text=fcall)))
      m<-nrow(Jn) # number of residuals / rows in Jacobian
      fcall<-paste(funjac,"(x,...)",sep='') # ?? ... handling
#      cat("fcall:",fcall,"\n")
      Ja<-try(eval(parse(text=fcall)))
      Jrd<-maxard10(Jn, Ja)
      if (cctrl$trace > 0) cat("Max rel diff shift 10 for Jacobian = ",Jrd,"\n")
      answer$Jrd<-Jrd        # If Jacobian second derivatives present, test code    
      funrsd<-paste(fname,".rsd",sep='')
      if(any(funrsd==ls(pattern=funrsd, name=".GlobalEnv"))) {
         fcall<-paste(funrsd,"(x,...)",sep='')
#         cat("fcall:",fcall,"\n")
         rsda<-try(eval(parse(text=fcall)))
         rsdn<-array(0, c(m,n,n)) # set up array for numerical Jacobian second derivatives
         # Create rsdn locally -- may want to adjust things
         hh<-1e-7 # Choose a fixed shift for derivs. Not necessarily a good idea??
         for (j in 1:n) {
            xadj<-x
            xadj[j]<-xadj[j]+hh
            fcall<-paste(funjac,"(xadj,...)",sep='')
            Jj<-try(eval(parse(text=fcall)))
            rsdn[,,j]<-(Jj-Ja)/hh
         }
         rsdrd<-maxard10(rsdn, rsda)
         if (cctrl$trace > 0) cat("Max rel diff shift 10 for Jacobian second derivatives = ",rsdrd,"\n")
         answer$rsdrd<-rsdrd
      } # end residual second derivative
    } # end jacobian stuff
  } # end of residuals
  # Test hessian if present
  funhess<-paste(fname,".h",sep='')
  answer$Hrd<-NA
  if(any(funhess==ls(pattern=funhess, name=".GlobalEnv"))) {
     fcall<-paste("hessian(",funf,",x,...)",sep='') # ?? ...
#     cat("Hessian fcall:",fcall,"\n")
     Hn<-try(eval(parse(text=fcall)))
#     cat("Hn:")
#     print(Hn)
#     cat("funhess:",funhess,"\n")
     fcall<-paste(funhess,"(x,...)",sep='')
     Ha<-try(eval(parse(text=fcall)))
     Hrd<-maxard10(Hn, Ha)
     if(cctrl$trace > 0) cat("Max rel diff shift 10 for Hessian = ",Hrd,"\n")
     answer$Hrd<-Hrd
  }
  return(answer) # ?? have not decided all content
}
### end of funcheck ***

