kktc <- function(par, fval, ngmax, nHes, 
     control=list(kkttol=1e-3, kkt2tol=1e-6, ktrace=FALSE)) {
# Provide a check on Kuhn-Karush-Tucker conditions based on quantities
# already computed. Some of these used only for reporting.
#
# Input:
#  par = a single vector of starting values
#  fval = objective function value
#  gvalmax = absolute value of the largest gradient component
#  hess = Hessian matrix evaluated at the parameters par
#
# Output: A list of four elements, namely,
#  ngmax = 
#  evratio = 
#  kkt1 = 
#  kkt2 = 
# 
# NOTE: Does not do projections for bounds and masks! 
#################################################################
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
#   print(control)
   if (is.null(control$ktrace)) trace<-FALSE else trace<-control$ktrace
   if (is.null(control$kkttol)) kkttol<-1e-3 else kkttol<-control$kkttol
   if (is.null(control$kkt2tol)) kkt2tol<-1e-6 else kkt2tol<-control$kkt2tol
   if (trace) { cat("kkttol=",kkttol,"   kkt2tol=",kkt2tol,"  trace=",trace,"\n") }
   if (trace) { cat("fval =",fval,"\n") }
   nbds<-0 # currently no bounds
   npar<-dim(nHes)[1]
   cat("Number of parameters =",npar,"\n")
   if (trace) { cat("KKT condition testing\n") }
   kkt1<-NA
   kkt2<-NA
   # test gradient
   if (trace) {
      cat("ngmax =",ngmax,"  test tol = ",kkttol*(1.0+abs(fval)),"\n")
   }
   kkt1<-(ngmax <= kkttol*(1.0+abs(fval)) ) # ?? Is this sensible?
   if (trace) {cat("KKT1 result = ",kkt1,"\n") }
# projections for bounds if they were present
#   bset<-sort(unique(c(which(par<=lower), which(par>=upper))))
#   nbds<-length(bset) # number of elements nulled by bounds
# Note: we assume that we are ON, not over boundary, but use <= and >=. No tolerance is used.
#   ng[bset]<-0 # "project" the gradient
#   nHes[bset,] <-0
#   nHes[, bset] <-0 # and the Hessian
#   ngatend <- ngatend
#   nhatend <- nhatend
#   cat("Hessian\n")
#   print(nHes)
   hev<- try(eigen(nHes)$values, silent=TRUE) # 091215 use try in case of trouble, # 20100711 silent
   if (trace) {
      cat("Hessian eigenvalues:\n")
      print(hev)
#      print("class(hev)=")
#      print(class(hev))
   }
   if (class(hev) != "try-error") {
# now look at Hessian
       negeig<-(hev[npar] <= (-1)*kkt2tol*(1.0+abs(fval))) # 20100711 kkt2tol
       evratio<-hev[npar-nbds]/hev[1]
       # If non-positive definite, then there will be zero eigenvalues (from the projection)
       # in the place of the "last" eigenvalue and we'll have singularity.
       # WARNING: Could have a weak minimum if semi-definite.
       kkt2<- (evratio > kkt2tol) && (! negeig) 
       if (trace) {cat("KKT2 result = ",kkt2,"\n") }
       cat("Hessian eigenvalues:\n")
       print(hev)
       cat("max abs gradient component = ",ngmax,"\n")
       cat("kkt1 =",kkt1,"    kkt2 =",kkt2,"\n")
       ans<-list(ngmax,evratio,kkt1,kkt2)
       names(ans)<-c("ngmax","evratio","kkt1","kkt2")
       return(ans)
   } else {
       stop("Eigenvalue failure")
   } # end kkt test
} ## end of optimx
