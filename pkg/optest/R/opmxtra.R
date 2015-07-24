opmxtra <- function(ansout, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            control=list(), ...) {

  cnames <- c(pstring, "value", "fevals", "gevals", "niter", "convergence", "kkt1", "kkt2", "xtimes")
  ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+8)

ans<-list(gmax,evratio,kkt1,kkt2,hev)
        names(ans)<-c("gmax","evratio","kkt1","kkt2","hev")

  cat("opm: wrapper to call optimr to run multiple optimizers\n")
   
  ans.ret <- list()

  allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf", 
                "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb")

  if (method == "ALL") control$all.methods <- TRUE
  if (control$all.methods) method <- allmeth # does not restrict if bounds??

  nmeth <- length(method)

# ?? do we need to restrict to bounded methods?? -- or leave to optimr??
  npar <- length(par)
  pstring<-names(par)
  if (is.null(pstring)) {
    pstring <- NULL
    for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
  } 
  cnames <- c(pstring, "value", "fevals", "gevals", "niter", "convergence", "kkt1", "kkt2", "xtimes")
  ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+8)
  colnames(ans.ret)<-cnames
  row.names(ans.ret)<-method
  cat("width of ans.ret =", npar+8,"\n")


  for (i in 1:nmeth) {
    meth <- method[i] # extract the method name
    # Note: not using try() here
    time <- system.time(ans <- optimr(par, fn, gr, method=meth, lower=lower, upper=upper, 
           hessian=hessian, control=control, ...))[1]
    cat("Method: ",meth,"\n")
    print(ans)
    # add to list
#    ans.ret <- c(ans.ret, c(method=meth, ans))

## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      if (control$trace>0) { cat("Post processing for method ",meth,"\n") }
      if (exists("ans$message")) {
           amsg<-ans$message
           ans$message <- NULL
      } else { amsg <- "none" }
      ngatend <- NA
      nhatend <- NA
      hev <- NA
      ans$gevals <- ans$counts[2]
      ans$fevals <- ans$counts[1]
      if ( control$save.failures || (ans$convergence < 1) ){# Save soln if converged or directed to save
#          if (control$trace && ans$convergence==0) cat("Successful convergence! \n") 
# Testing final soln. Use numDeriv for gradient & Hessian; compute Hessian eigenvalues
          hessOK<-FALSE # indicator for later
          gradOK<-FALSE
#          if ((control$kkt || hessian) && (ans$convergence != 9999)) {
#              if (control$trace>0) cat("Compute Hessian approximation at finish of ",method[i],"\n")
#              if (!is.null(uhess)){ # check if we have analytic hessian 
#                 nhatend<-try(uhess(ans$par, ...), silent=TRUE)
#                 if (class(nhatend) != "try-error") {
#                    hessOK<-TRUE
#                 }
#              } else {
#                 if (is.null(ugr)) {
#                     nhatend<-try(hessian(ufn, ans$par, ...), silent=TRUE) # change 20100711
#                 } else {
#                     nhatend<-try(jacobian(ugr,ans$par, ...), silent=TRUE) # change 20100711
#                 } # numerical hessian at "solution"
#                 if (class(nhatend) != "try-error") {
#                    hessOK<-TRUE
#                 }
#              } # end hessian calculation
#          } # end test if hessian computed
          ans$kkt1<-NA
          ans$kkt2<-NA
#          if ((hessian || control$kkt) && (ans$convergence != 9999)) {# avoid test when method failed
#             if (control$trace>0) cat("Compute gradient approximation at finish of ",method[i],"\n")
#             if (is.null(ugr)) {
#                 ngatend<-try(grad(ufn, ans$par, ...), silent=TRUE) # change 20100711
#             } else {
#                 ngatend<-try(ugr(ans$par, ...), silent=TRUE) # Gradient at solution # change 20100711
#             }
#             if (class(ngatend) != "try-error") gradOK<-TRUE # 100215 had == rather than != here
#             if ( (! gradOK) && (control$trace>0)) cat("Gradient computation failure!\n") 
#             if (gradOK) {
#                # test gradient
#                ans$kkt1<-(max(abs(ngatend)) <= control$kkttol*(1.0+abs(ans$value)) ) # ?? sensible?
#                if (hessOK) {
#                   # For bounds constraints, we need to "project" the gradient and Hessian
#                   bset<-sort(unique(c(which(ans$par<=lower), which(ans$par>=upper))))
#                   nbds<-length(bset) # number of elements nulled by bounds
#                   # Note: we assume that we are ON, not over boundary, 
#                   # but use <= and >=. No tolerance is used.
#                   ngatend[bset]<-0 # "project" the gradient
#                   nhatend[bset,] <-0
#                   nhatend[, bset] <-0 # and the Hessian
#                   if (!isSymmetric(nhatend, tol=sqrt(.Machine$double.eps))) {
#                      # hessOK<-FALSE ## Assume we will keep hessian after symmetrizing
#                      asym <- sum(abs(t(nhatend) - nhatend))/sum(abs(nhatend))
#                      asw <- paste("Hessian is reported non-symmetric with asymmetry ratio ", 
#                      asym, sep = "")
#                      if (control$trace > 1) cat(asw, "\n")
#                      if (control$dowarn) warning(asw)
#                      ### if (asym > control$asymtol) stop("Hessian too asymmetric") ##??as yet don't stop
#                      if (control$trace > 1) cat("Force Hessian symmetric\n")
#                      if (control$dowarn) warning("Hessian forced symmetric", call. = FALSE)
#                      nhatend <- 0.5 * (t(nhatend) + nhatend)
#                   }  # end symmetry test
#                   hev<- try(eigen(nhatend)$values, silent=TRUE) # 091215 use try in case of trouble
#                   if (control$kkt){
#   	              if (class(hev) != "try-error") {# answers are OK, check Hessian properties
#                         if (any(is.complex(hev))){
#                            hessOK<-FALSE
#                            cat("Complex eigenvalues found for method =",meth,"\n")
#                            cat("coefficients for function value", ans$value," :\n")
#                            print(ans$par)
#                            dput(nhatend, file="badhess.txt")
#                            warning("Complex eigenvalues found for method =",meth)
#                         }
#                         if (hessOK) {
#                            negeig<-(hev[npar] <= (-1)*control$kkt2tol*(1.0+abs(ans$value))) 
#                            evratio<-hev[npar-nbds]/hev[1]
#                            # If non-positive definite, then there have zero evs (from the projection)
#                            # in the place of the "last" eigenvalue and we'll have singularity.
#                            # WARNING: Could have a weak minimum if semi-definite.
#                            ans$kkt2<- (evratio > control$kkt2tol) && (! negeig)
#                         }
#                      } else {
#                         warnstr<-paste("Eigenvalue failure after method ",method[i],sep='')
#                         if (control$dowarn) warning(warnstr)
#                         if (control$trace>0) {
#                            cat("Hessian eigenvalue calculation failure!\n")
#                            print(nhatend)
#                         }
#                      }
#                   } # kkt2 evaluation
#                } else { # computing Hessian has failed
#                   warnstr<-paste("Hessian not computable after method ",method[i],sep='')
#                   if (control$dowarn) warning(warnstr)
#                   if (control$trace>0) cat(warnstr,"\n") 
#                }
#             } else { # gradient failure
#                warnstr<-paste("Gradient not computable after method ",method[i],sep='')
#                if (control$dowarn) warning(warnstr)
#                if (control$trace>0) cat(warnstr,"\n") 
#             }
#          } # end kkt test
          cat("str(time):\n")
          print(time)
          ans$xtimes <- time
          # Do we want more information saved?
          if (control$trace>0) { 
		cat("Save results from method ",meth,"\n") 
	  	print(ans)
	  }
	  if (control$trace>0) { cat("Assemble the answers\n") }
          cat("ans.ret now\n")
          print(ans.ret)
          ans$nitns <- NA
          cat("add ans for meth = ",meth,"\n")
          print(ans)
          addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, ans$nitns,
                              ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
          cat("length addvec = ",length(addvec),"\n")
          ans.ret[meth, ] <- addvec
          if (! gradOK) ngatend <- NA
          if (! hessOK) {
             nhatend <- NA
             hev <- NA
          }
      }  ## end post-processing of successful solution
#      ans.details<-rbind(ans.details, list(method=meth, ngatend=ngatend, nhatend=nhatend, hev=hev, message=amsg))
      # 1303234 try list() not c()
#      row.names(ans.details)[[i]]>=meth
#      	if (control$follow.on) {
#		par <- ans$par # save parameters for next method
#		if (i < nmeth && (control$trace>0)) cat("FOLLOW ON!\n") # NOT trace ??
#	}
    } # End loop  ## end loop over method (index i)
    ansout <- NULL # default if no answers
    if (length(ans$par) > 0) { # cannot save if no answers
	ansout <- data.frame(ans.ret)# Don't seem to need drop=FALSE
#        attr(ansout, "details")<-ans.details
    }
    ansout # return(ansout)

} ## end of optimx

