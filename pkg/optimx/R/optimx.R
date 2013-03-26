optimx <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

  optcfg <- optimx.setup(par, fn, gr, hess, lower, upper, 
            method, itnmax, hessian, control, ...)
# Parse and use optcfg
  if (optcfg$ctrl$starttests) {
    optchk <- optimx.check(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, lower,
           upper, hessian, optcfg$ctrl, have.bounds=optcfg$have.bounds, ...)
  }
  ansout <- optimx.run(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, lower, upper,
            optcfg$method, itnmax, hessian, optcfg$ctrl, ...)
  details <- attr(ansout, "details")
  attr(ansout, "details") <- NULL ## Not sure this necessary
  if (optcfg$ctrl$maximize) {
     if (control$trace) cat("Reversing sign on objective, gradient, & hessian\n")
     ansout$value <- - ansout$value
     nlist<-dim(details)[[1]]
     for (i in 1:nlist) {
        details[[i,"ngatend"]] <- - details[[i,"ngatend"]] 
        details[[i,"nhatend"]] <- - details[[i,"nhatend"]] 
        details[[i,"hev"]] <- - details[[i,"hev"]] 
     }
  }
  structure(ansout, details = details, maximize = optcfg$ctrl$maximize,
            npar = optcfg$npar, class = c("optimx", "data.frame"))
} ## end of optimx

