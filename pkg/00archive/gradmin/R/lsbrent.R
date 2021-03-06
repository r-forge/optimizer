lnbrent<-function(fn, fbest, xc, d, grv, w, ...) { # Line search using internal optimize()
  cat("fn:\n")
  print(fn)
  ## Uses Brent's method to find the best stepsize in interval
  # fbest is best function value so far. NOT used.
  # grv is numeric gradient vector -- NOT used
  # ?? more documentation
  flsch<-function(st) {
    # computes the function value at stepsize st on line (xc + gm*d)
    # Essentially flsch(st)
    # gm: step size
    # fn: objective function
    # xc: base set of parameters
    # d : search direction
    #      nf <- nf +1
    fval<-fn(xc+st*d,...)
    fval
  }
  cat("function at ends of interval\n")
  sta <- w$stepmin
  cat("f(",sta,")=", flsch(sta),"\n")
  stb <- w$stepmax
  cat("f(",stb,")=", flsch(stb),"\n")
  
  #  lout<-optimize(flsch,interval=c(w$stepmin, w$stepmax),
  #                  lower=w$stepmin, upper=w$stepmax,...)
  # note fmin rather than objective in return  
  lout<-pracma:fminbnd(flsch,w$stepmin, w$stepmax, ...)
  cat("lnsrch lout:")
  print(lout)
  rlout <- lout$xmin
  #  cat("structure of rlout")
  #  print(str(rlout))
  attr(rlout, "Fval") <- lout$fmin
  attr(rlout, "fcount") <- (lout$niter + 1) # fevals is iterations + 1
  rlout # Note: returns stepsize, not x
} # end default line search
