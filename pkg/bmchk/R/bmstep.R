bmstep <- function(par, srchdirn, lower=NULL, upper=NULL, bdmsk=NULL, trace=0) {  
## Find maximum step from par along srchdirn given bounds and masks
# 20101231 ?? check use of par and bvec -- can we simplify?

# 20101031 -- issue of bounds not working correctly
#  - -Inf seems to upset bounds
#  - no upper bounds gives troubles (applies to Rcgmin too!)
#
# Input:
#  par = a vector containing the starting point
#  srchdirn = the direction of search (a vector)
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
#    A double giving the maximum step. Not bigger than maxstep.
#
########## length of vectors #########
n<-length(par)
############# bounds and masks ################
# check if there are bounds
  if(is.null(lower) || ! any(is.finite(lower))) nolower<-TRUE else nolower<-FALSE
  if(is.null(upper) || ! any(is.finite(upper))) noupper<-TRUE else noupper<-FALSE
# Next line NOT same as in bmchk(). Leave out bdmsk.
  if(nolower && noupper) bounds<-FALSE else bounds<-TRUE
  if (any(bdmsk==0)) {
     if (trace > 2) cat("Masks present -- adjusting search direction.\n")
     srchdirn[which(bdmsk==0)]<-0 # adjust search direction for masked elements
  }
  if(nolower) lower<-rep(-Inf,n)
  if(noupper) upper<-rep(Inf,n)
######## find maximum step (may be Inf) #############
# distance to bounds
  # test bmstep ideas
  sn<-rep(Inf, n)
  sp<-sn
# these are the max steps to bounds
  idxn<-which(srchdirn<0)
  idxp<-which(srchdirn>0)
# try out these
  suppressWarnings(sn[idxn]<-(par-lower)/srchdirn)
  suppressWarnings(sp[idxp]<-(upper-par)/srchdirn)
  maxstep<-min(sn,sp)
} ## end of bmstep.R
