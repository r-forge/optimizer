opm <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

  cat("opm: wrapper to call optimr to run multiple optimizers\n")
  npar<-length(par)
  if (is.null(control$trace)) control$trace <- 0
  if (control$trace > 1) cat("Function has ",npar," arguments\n")

# ?? check before or after scaling -- better before, but may want sometime
#    to be able to check my transformations

  if (is.null(control$starttests) || control$starttests) {
     # Check parameters are in right form
     if (!is.null(dim(par))) stop("Parameters should be a vector, not a matrix!")
     if (! is.vector(par) ) { stop("The parameters are NOT in a vector")  }


     # Check parameters in bounds (090601: As yet not dealing with masks ??)
     infeasible<-FALSE
#     Even if we have control$have.bounds, do the check
     if (is.null(control$keepinputpar)) { shift2bound <- TRUE } 
     else {shift2bound <- ! control$keepinputpar}
     bc <- bmchk(par, lower=lower, upper=upper, trace=control$trace, shift2bound)
     if (! bc$admissible) stop("At least one lower bound is > corresponding upper bound")
     if (infeasible && control$dowarn) warning("Parameters may be out of bounds")
     if (control$trace > 0) {
        cat("Parameter relation to bounds\n")
        print(bc$bchar)
     }
     if (bc$parchanged) par <- bc$bvec # adjust parameters to bounds
     control$have.bounds <- bc$bounds # reset the control vector element for bounds

     # Check if function can be computed
     tmp <- readline("about to call fnchk")
     cat("Call fnchk with par:")
     print(par)
     cat("Call fnchk with fn:")
     print(fn)
     checkfn <- fnchk(par, fn, trace=control$trace, ...)

     if (checkfn$infeasible) {
        cat("fnchk exit code and msg:",checkfn$excode," ",checkfn$msg,"\n")
        stop("Cannot evaluate function at initial parameters")
     }

     grbad <- FALSE
     if (is.null(control$usenumDeriv)) control$usenumDeriv <- FALSE # ?? may be not correct
     if (! is.null(gr) && ! is.character(gr) && ! control$usenumDeriv){ # check gradient
     #   have gr fn        not specified method     not using numDeriv
       gname <- deparse(substitute(gr))
       if (control$trace > 0) cat("Analytic gradient from function ",gname,"\n\n")
       fval <- ufn(par, ...) 
       gn <- grad(func=fn, x=par,...) # 
       ga <- gr(par, ...)
       # Now test for equality (090612: ?? There may be better choices for the tolerances.
       teps <- (.Machine$double.eps)^(1/3)
       grtest <- max(abs(gn-ga))/(1 + abs(fval))
       if (control$trace > 0) cat("Gradient test value =",grtest,"\n")
       if (grtest >= teps) {
         stop("Gradient function test value ",grtest," check gradient! \n", call.=FALSE)
         grbad <- TRUE # Never get here if we stop ??
       }
       } else if (control$trace > 0) cat("Analytic gradient not made available.\n")

       hessbad <- FALSE
##?? need to fix this
       if (! is.null(hess) && ! is.character(hess)){ # check Hessian - if character then numeric
          hname <- deparse(substitute(uhess))
          if (control$trace > 0) cat("Analytic hessian from function ",hname,"\n\n")
          hn <- hessian(func=ufn, x=par, parscale=control$parscale, ...) # ?? should we use dotdat
          ha <- uhess(par, ...)
          # Now test for equality
          teps <- (.Machine$double.eps)^(1/3)
          if (max(abs(hn-ha))/(1 + abs(fval)) >= teps) stop("Hessian function might be wrong - check it! \n", call.=FALSE)
          optchk$hessbad <- TRUE
       } else if (control$trace > 0) cat("Analytic Hessian not made available.\n")

       # Scaling check  091219
       cat("scale check\n")
       scalebad <- FALSE
       srat<-scalecheck(par, lower, upper,control$dowarn)
       sratv<-c(srat$lpratio, srat$lbratio)
       if (is.null(control$scaletol)) control$scaletol <- 3
       if (max(sratv,na.rm=TRUE) > control$scaletol) { 
	  warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
          if (! is.null(control$dowarn) && control$dowarn)  warning(warnstr)
       }
       if (control$trace > 0) {
         cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
       }

  
  } # end starttests

tmp <- readline("Run optest.setup")

  optcfg <- optest.setup(par, fn=fn, gr=gr, hess=hess,
            lower=lower, upper=upper, method=method, 
            itnmax=itnmax, hessian=hessian, control=control, ...)
  print(str(optcfg))
tmp <- readline("Now run optimr")

  tctrl <- optcfg$ctrl
  tctrl$parscale <- NULL # make sure we don't have scaling double called
  # the scaling is passed via pscale for ufn, ugr and (possibly uhess??)
  print(str(tctrl))

#  if (optcfg$ctrl$have.bounds) {
  ans <- optimr(par=optcfg$spar, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, 
         hessian=hessian, control=tctrl, pscale=control$parscale, ...)
#  }    
  ans$par <- ans$par * control$parscale
  ans

} ## end of optimx

