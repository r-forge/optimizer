optest <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...) {

#  cat("optest: control is \n")
#  print(control)
#  tmp <- readline("cont. optest")

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
#  opchk <- optest.check(par, ufn=optcfg$ufn, ugr=optcfg$ugr, lower=lower, upper=upper, 

  opchk <- optest.check(par, ufn=optcfg$ufn, ugr=optcfg$ugr, 
         ctrl=optcfg$ctrl, ...)
#  print(opchk)
#  tmp <- readline("try optimr")
#  cat("optimr call for method ",optcfg$method," using trace =", optcfg$ctrl$trace,"\n")

  if (optcfg$ctrl$have.bounds) {
  ans <- optimr(par, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, lower=lower, upper=upper, 
         hessian=hessian, control=optcfg$ctrl, ...)
  } else {
  ans <- optimr(par, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, 
         hessian=hessian, control=optcfg$ctrl, ...)
  }    
#  print(ans)
#  tmp <- readline("try optest.run")
#  ans <- optim(par, fn=optcfg$ufn, gr=optcfg$ugr, method=optcfg$method, lower=lower, upper=upper, 
#         hessian=hessian, control=optcfg$ctrl, ...)
#  print(ans)
  ans    



} ## end of optimx

