bmchk <- function(par, lower=NULL, upper=NULL, bdmsk=NULL, trace=0) {  
## Bounds and masks check
#  ?? check use of par and bvec -- can we simplify?

# 20101031 -- issue of bounds not working correctly
#  - -Inf seems to upset bounds
#  - no upper bounds gives troubles (applies to Rcgmin too!)
#
# Input:
#  par = a vector containing the starting point
#  lower = vector of lower bounds on parameters
#  upper = vector of upper bounds on parameters
#    Note: free parameters outside bounds will be adjusted to bounds.
#  bdmsk = control vector for bounds and masks. Parameters for which bdmsk are 1
#         are unconstrained or "free", those with bdmsk 0 are masked i.e., fixed.
#         For historical reasons, we use the same array as an indicator that a
#         parameter is at a lower bound (-3) or upper bound (-1)
#  trace = control of output: 0 for none (default), >0 for output
##
# Output:
#    A list with components: 
#     bvec: The best set of parameters found.
#     bdmsk: adjusted inputs
#     nolower: TRUE if no lower bounds, FALSE otherwise
#     noupper: TRUE if no upper bounds, FALSE otherwise
#     bounds:  TRUE if any bounds, FALSE otherwise
#     admissible: TRUE if admissible, FALSE if not
#     parchanged: TRUE if parameters changed, FALSE if not
#
############# bounds and masks ################
# set default masks if not defined
   n<-length(par)
   bvec<-par # copy parameters (will this impact memory use??)
   if (is.null(bdmsk)) {
      bdmsk<-rep(1,n)
   }
   if (trace > 2) {
      cat("bdmsk:")
      print(bdmsk)
   }
# check if there are bounds
  if(is.null(lower) || ! any(is.finite(lower))) nolower<-TRUE else nolower<-FALSE
  if(is.null(upper) || ! any(is.finite(upper))) noupper<-TRUE else noupper<-FALSE
  if(nolower && noupper && all(bdmsk == 1)) bounds<-FALSE else bounds<-TRUE
  if (trace > 2) cat("Bounds: nolower = ",nolower,"  noupper = ",noupper," bounds = ",bounds,"\n")
  if(nolower) lower<-rep(-Inf,n)
  if(noupper) upper<-rep(Inf,n)

######## check bounds and masks #############
## NOTE: do this inline to avoid call (??should we change this?)
  if (bounds) {
     # Make sure to expand lower and upper
     if(! nolower & (length(lower)<n)) {
        if (length(lower)==1) { lower<-rep(lower,n) } else { stop("1<length(lower)<n") }
     } # else lower OK
     if(! noupper & (length(upper)<n)) {
        ## tmp<-readline("Check length upper ")
        if (length(upper)==1) { upper<-rep(upper,n) } else { stop("1<length(upper)<n") }
     } # else upper OK
     # At this point, we have full bounds in play
     ######## check admissibility ########
     admissible<-TRUE
     parchanged<-FALSE
     if (any(lower[which(bdmsk != 0)] > upper[which(bdmsk != 0)])) admissible <- FALSE
     tol<-.Machine$double.eps*max(abs(upper),abs(lower),1)
     maskadded<-FALSE
     if (any((upper-lower)<tol)) { # essentially masked
           warning("One or more lower bounds equals an upper bound, imposing mask")
           bdmsk[which(upper-lower < tol)] = 0
           maskadded<-TRUE
     }
     if (trace > 0) cat("admissible = ",admissible,"\n")
     if (admissible) {
        # This implementation as a loop, but try later to vectorize
        for (i in 1:n) {
           if (bdmsk[i] == 0) { 
              # NOTE: we do not change masked parameters, even if out of bounds
              if(! nolower) {
                 if((bvec[i]<lower[i]) & (trace > 0)) { 
                    cat("WARNING: ",bvec[i]," = MASKED x[",i,"] < lower bound = ",
                        lower[i],"\n")
                 }
              }
              if(! noupper) {
                 if ((bvec[i]>upper[i]) & (trace > 0)) { 
                    cat("WARNING: ",bvec[i]," = MASKED x[",i,"] > upper bound = ",
                        upper[i],"\n")
                 }
              }
           } else { # not masked, so must be free or active constraint
              if(! nolower){
                 if (bvec[i]<=lower[i]) { # changed 090814 to ensure bdmsk is set
		   if (trace > 0) {
                      cat("WARNING: x[",i,"], set ",bvec[i]," to lower bound = ",lower[i],"\n")
                   }
                   bvec[i]<-lower[i]
                   parchanged<-TRUE
                   bdmsk[i] <- -3 # active lower bound
                 }
              }
              if(! noupper) {
                 if (bvec[i]>=upper[i]) { # changed 090814 to ensure bdmsk is set 
		   if (trace > 0) {
                      cat("WARNING: x[",i,"], set ",bvec[i]," to upper bound = ",upper[i],"\n")
                   }
                   bvec[i]<-upper[i]
                   parchanged<-TRUE
                   bdmsk[i] <- -1 # active upper bound
                 }
              }
           } # end not masked
        } # end loop for bound/mask check
     }
   }
   ############## end bounds check #############
   bcout<-list(bvec, bdmsk, nolower, noupper, bounds, admissible, maskadded, parchanged)
   names(bcout)<-c("bvec", "bdmsk", "nolower", "noupper", "bounds", "admissible", 
        "maskadded", "parchanged")
   return(bcout)
} ## end of bmchk.R
