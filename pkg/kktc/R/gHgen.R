gHgen <- function(par, fn, gr=NULL, hess=NULL, ktrace=FALSE, 
            control=list(asymtol=1.0e-7), ...) {
#  Generate the gradient and Hessian for a given function at the parameters par.
#
# Input:
#  par = a vector for function parameters
#  fn = a user function (assumed to be sufficeintly differentiable)
#  gr = name of a function to compute the (analytic) gradient of the user function
#  hess = name of a function to compute the (analytic) hessian of the user function
#         This will rarely be available, but is included for completeness.
#  ktrace = logical flag (default FALSE) to monitor progress
#  control=list of controls. Currently only asymtol=1.0e-7.
#         This is used to decide if we should force symmetry.
#  ...    = additional arguments to the objective, gradient and Hessian functions
#
# Output:
# A list containing
#    g = a vector containing the gradient
#    H = a matrix containing the Hessian
#
#  Author: John Nash
#  Date:  January 14, 2011
#
#################################################################
   asymtol<-control$asymtol # used to judge if we should symmetrize Hessian
   require(numDeriv)
   gradOK<-FALSE
   if (ktrace) cat("Compute gradient approximation\n")
   if (is.null(gr)) {
      ng<-try(grad(fn, par, ...), silent=TRUE) # change 20100711
   } else {
      ng<-try(gr(par, ...), silent=TRUE) # Gradient at solution # change 20100711
   }
   if (class(ng) != "try-error") gradOK<-TRUE # 100215 had == rather than != here
   if ( ! gradOK ) stop("Gradient computations failure!")
   if (ktrace) print(ng)
   if (ktrace) cat("Compute Hessian approximation\n")
   if (is.null(hess)) {
      if (is.null(gr)) {
         nH<-try(hessian(fn, par, ...), silent=TRUE) # change 20100711
         if (class(nH) == "try-error") stop("Unable to compute Hessian using numDeriv::hessian")
	 # Do not need to check for symmetry
      } else {
         nH<-try(jacobian(gr,par, ...), silent=TRUE) # change 20100711
         if (class(nH) == "try-error") stop("Unable to compute Hessian using numderiv::jacobian")
         if (! isSymmetric(nH) ) {
            asym<-sum(abs(t(nH)-nH))/sum(abs(nH))
            asw<-paste("nH from jacobian is reported non-symmetric with asymmetry ratio ",asym,sep='')
            if (ktrace) cat(asw,"\n")
            warning(asw)
            if (asym > asymtol) stop("Hessian too asymmetric")
            if (ktrace) cat("Force Hessian symmetric\n")
            else warning("Hessian forced symmetric")
            nH<-0.5*(t(nH)+nH)
         } # end if ! isSymmetric
      } # numerical hessian at "solution"
   } else { 
      nH<-try(hess(par,...), silent=TRUE) # change 20110222
      if (class(nH)=="try-error") stop("Hessian evaluation with function hess() failed")
      if (! isSymmetric(nH) ) {
         asym<-sum(abs(t(nH)-nH))/sum(abs(nH))
         asw<-paste("nH from hess() is reported non-symmetric with asymmetry ratio ",asym,sep='')
         if (ktrace) cat(asw,"\n")
         warning(asw)
         if (asym > asymtol) stop("Hessian too asymmetric")
         if (ktrace) cat("Force Hessian symmetric\n")
         else warning("Hessian forced symmetric")
         nH<-0.5*(t(nH)+nH)
      } # end if ! isSymmetric
   } # end hessian computation
   if (ktrace) print(nH)
   ansout<-list(ng, nH)
   names(ansout)=c("ng","nH")
   return(ansout)
} ## end gHgen
