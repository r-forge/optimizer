rootsearch <- function(fn=fn, start=NULL, ftrace=TRUE, init.step=1, stepmult = -1.5,
                 llmin=-1e+7, uumax=1e+7, ...){
  ## A routine to search for an interval which brackets a root for 
  ## function fn(x, ...) 
  # fn = the function for which root is to be found
  # start = value of x to use as first try
  # ftrace = logical TRUE if method is to be traced
  # init.step = value to add to init.start for next try
  # stepmult = multiplier to adjust step
  # llmin, uumax: extreme limits. Stop search if we exceed either.
  # ... = additional arguments to fn or FnTrace (e.g., ifn)
  
  ## OUTPUT
  #  llow: lower end of interval
  #  uup : upper end of interval
  #  flow: function value at llow
  #  fup : function value at uup
  #  count: number of function evaluations used
  #         Set negative if we have failed.

##  dots <- list(...)
  if (init.step == 0) stop("Cannot have zero step")
  kf <- 2
  x <- start
  step <- init.step
  if (step < 0) {
     x <- start + init.step
     step <- -step 
  } # end reset for negative step
  llow <- x 
  flow <- fn(llow, ...)
  uup <- x + step
  x <- uup
  fup <- fn(uup, ...)
  while (sign(flow)*sign(fup) >= 0) { ## Note sign 0 can be 0
     if(ftrace)  cat("f(",llow,")=",flow,"  f(",uup,")=",fup,"   step=",step,"\n")
      step <- stepmult * step # increase size and change sign
      x <- x + step
      if ((x < llmin) || (x > uumax) ) {
          kf <- - kf
          break # failed
      } 
      kf <- kf + 1
      if (x < llow) {
         llow <- x
         flow <- fn(llow, ...)
      } 
      else {
         uup <- x
         fup <- fn(uup, ...)
      }
  } # end while
  res <- list(llow=llow, uup = uup, flow = flow, fup = fup, count = kf)
}


f1 <- function(x) { 2*(x - 10) }

cat("f1 try\n")
res1 <- rootsearch(f1, 1)
print(res1)

f2 <- function(x) (0.5*(x-1)^2)
res2 <- rootsearch(f2, 1)
print(res2)


