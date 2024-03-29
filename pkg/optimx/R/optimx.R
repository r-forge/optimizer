optimx <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {
## NOTE: optimx is NOT easily maintainable, but is here for legacy purposes.
##  Please use opm, optimr, polyopt and multistart.
  if (is.null(control$trace)) control$trace<-0
  optcfg <- optimx.setup(par, fn, gr, hess, lower, upper, 
            method, itnmax, hessian, control, ...)
# Parse and use optcfg
  if (optcfg$ctrl$starttests) {
    if (! is.null(gr)) {
       tg<-grchk(par, ffn=fn, ggr=gr, trace=control$trace, ...)
       if (! tg) stop("optimx: Gradient function may be incorrect.")
    }
    if (! is.null(hess)) {
       th<-hesschk(par, ffn=fn, ggr=gr, hhess=hess, trace=control$trace, ...)
       if (! th) stop("optimx: Hessian function may be incorrect.")   
    }
    srat<-scalechk(par, lower, upper, dowarn=TRUE)
    print(srat) # Always print if starttests selected. Overrides trace.
    sratv<-c(srat$lpratio, srat$lbratio)
    if (max(sratv,na.rm=TRUE) > 3) { # scaletol from ctrldefault in optimx
      warnstr<-"Parameters or bounds appear to have different scalings.\n
      This can cause poor performance in optimization. \n
      It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
      cat(warnstr,"\n")
    }
  }
  optcfg$ctrl$have.bounds<-optcfg$have.bounds # to pass boundedness
  if (! is.null(control$trace) && control$trace > 1) {
    cat("optcfg:"); print(optcfg)
  }

  ansout <- optimx.run(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, lower, upper,
            optcfg$method, itnmax, hessian, optcfg$ctrl)
  details <- attr(ansout, "details")
  attr(ansout, "details") <- NULL ## Not sure this necessary, but null here and replace below
  if (optcfg$ctrl$maximize) {
     if (optcfg$ctrl$trace>0) cat("Reversing sign on objective, gradient, & hessian\n")
     ansout$value <- - ansout$value
     nlist<-dim(details)[[1]]
     for (i in 1:nlist) {
        details[[i,"ngatend"]] <- - details[[i,"ngatend"]] 
        details[[i,"nhatend"]] <- - details[[i,"nhatend"]] 
        details[[i,"hev"]] <- - details[[i,"hev"]] 
     }
  }
  rownames(details) <- details[, "method"]
  ##JN -- don't remove method: 
  ##JN  details <- details[, colnames(details) != "method", drop=FALSE]
  # Fix kkt test output to logical
    ansout[ , "kkt1"] <- as.logical(ansout[ , "kkt1"])
    ansout[ , "kkt2"] <- as.logical(ansout[ , "kkt2"])

  answer <- structure(ansout, details = details, maximize = optcfg$ctrl$maximize,
            npar = optcfg$npar, follow.on=optcfg$ctrl$follow.on,
            class = c("optimx", "data.frame"))
  answer # requested by Gabor 1408
} ## end of optimx

