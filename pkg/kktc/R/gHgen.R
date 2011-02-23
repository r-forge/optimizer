gHgen <- function(par, fn, gr=NULL, hess=NULL, ktrace=FALSE, ...) {
#  Generate the gradient and Hessian for a given function at the parameters par.
#
# Input:
#  par = a vector for function parameters
#  fn = a user function (assumed to be sufficeintly differentiable)
#  gr = name of a function to compute the (analytic) gradient of the user function
#  hess = name of a function to compute the (analytic) hessian of the user function
#         This will rarely be available, but is included for completeness.
#  ktrace = logical flag (default FALSE) to monitor progress
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
      } else {
         nH<-try(jacobian(gr,par, ...), silent=TRUE) # change 20100711
      } # numerical hessian at "solution"
   } else { 
         nH<-try(hess(par,...), silent=TRUE) # change 20110222
   } # end hessian computation
   if (class(nH) == "try-error") stop("Unable to compute Hessian")
   if (ktrace) print(nH)
   ansout<-list(ng, nH)
   names(ansout)=c("ng","nH")
   return(ansout)
} ## end gHgen


