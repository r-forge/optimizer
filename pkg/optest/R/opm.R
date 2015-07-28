opm <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=FALSE,
            control=list(),
             ...) {

  npar <- length(par)
  pstring<-names(par)
  ctrl <- ctrldefault(npar)
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
  }
  control <- ctrl
  if(control$trace > 0) cat("opm: wrapper to call optimr to run multiple optimizers\n")

  allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf", 
                "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb")

  if (length(method) == 1 && method == "ALL") control$all.methods <- TRUE
  if (control$all.methods) method <- allmeth # does not restrict if bounds??

  nmeth <- length(method)

# ?? do we need to restrict to bounded methods?? -- or leave to optimr??
  if (is.null(pstring)) {
    pstring <- NULL
    for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
  } 
#  cnames <- c(pstring, "value", "fevals", "gevals", "niter", "convergence", "kkt1", "kkt2", "xtimes")
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence", "kkt1", "kkt2", "xtimes")
  ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+7)
  ans.ret <- data.frame(ans.ret)
  colnames(ans.ret)<-cnames
  row.names(ans.ret)<-method
  ans.details <- list()
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
          if ((control$trace > 0) && (ans$convergence==0)) cat("Successful convergence! \n") 
# Testing final soln. Use numDeriv for gradient & Hessian; compute Hessian eigenvalues
          if ((control$kkt || hessian) && (ans$convergence != 9999)) {
             kktres <- kktchk(ans$par, fn, gr, hess=NULL, upper=NULL, lower=NULL, 
                    maxfn=control$maximize, control=control, ...) 
             ans$kkt1<-as.logical(kktres$kkt1)
             ans$kkt2<-as.logical(kktres$kkt2)
          }
# put together results
#          cat("str(time):\n")
#          print(time)
          ans$xtimes <- time
          # Do we want more information saved?
          if (control$trace > 0) { 
		cat("Save results from method ",meth,"\n") 
	  	print(ans)
	  }
	  if (control$trace>0) { cat("Assemble the answers\n") }
#          cat("ans.ret now\n")
#          print(ans.ret)
#          ans$nitns <- NA
#          cat("add ans for meth = ",meth,"\n")
#          print(ans)
          addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, 
                              ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
#          cat("length addvec = ",length(addvec),"\n")
          ans.ret[meth, ] <- addvec
      }  ## end post-processing of successful solution
      ans.details<-rbind(ans.details, list(method=meth, ngatend=kktres$ngatend, 
             nhatend=kktres$nhatend, hev=kktres$hev, message=amsg))
      # 1303234 try list() not c()
      row.names(ans.details)[[i]]<-meth
    } # End loop  ## end loop over method (index i)
    ansout <- NULL # default if no answers
    if (length(ans$par) > 0) { # cannot save if no answers
	ansout <- ans.ret # Don't seem to need drop=FALSE
        attr(ansout, "details")<-ans.details
        ansout[, "kkt1"] <- as.logical(ansout[, "kkt1"])
        ansout[, "kkt2"] <- as.logical(ansout[, "kkt2"])
    }
    ansout # return(ansout)
} ## end of opm

