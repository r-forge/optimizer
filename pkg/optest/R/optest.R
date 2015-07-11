optest <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

#  cat("optest: control is \n")
#  print(control)
#  tmp <- readline("cont. optest")

  ansu <- optimr(par, fn, gr, method=method, control=control, ...) # unscaled

  optcfg <- optest.setup(par, fn=fn, gr=gr, hess=hess,
            lower=lower, upper=upper, method=method, 
            itnmax=itnmax, hessian=hessian, control=control, ...)
#  print(optcfg)

#  cat("Test the optcfg function and gradient:\n")
#  print(optcfg$ufn)
#  print(optcfg$ugr)
#  fval <- optcfg$ufn(par, ...)
#  gval <- optcfg$ugr(par, ...)
#  cat("fval at start is ",fval,"\n")
#  cat("gval at start is ")
#  print(gval)

#  tmp <- readline("try optest.check")
# ?? need to add bounds here for checking!!! Note bounds will be scaled too!!!
# ?? temp shut of check
#  opchk <- optest.check(par=optcfg$spar, ufn=optcfg$ufn, ugr=optcfg$ugr, ctrl=optcfg$ctrl, ...)
#  print(opchk)
#  tmp <- readline("try optimr")
#  cat("optimr call for method ",optcfg$method," using trace =", optcfg$ctrl$trace,"\n")

  tctrl <- optcfg$ctrl
  tctrl$parscale <- NULL # make sure we don't have scaling


#  if (optcfg$ctrl$have.bounds) {
#  ans <- optimr(par=optcfg$spar, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, 
#	lower=lower, upper=upper, hessian=hessian, control=tctrl,  pscale=control$parscale, ...)
#  } else {
  ans <- optimr(par=optcfg$spar, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, 
         hessian=hessian, control=tctrl, pscale=control$parscale, ...)
#  }    
  list(ans=ans, ansu=ansu)    # ??? still need to unscale

} ## end of optimx

