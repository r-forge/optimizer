# hjn.R -- R implementation of J Nash BASIC HJG.BAS 20160705
hjn <- function(par, fn, lower=-Inf, upper=Inf, bdmsk=NULL, control=list(trace=0), ...){
  n <- length(par) # number of parameters
  if (is.null(control$trace)) control$trace <- 0 # just in case
  if (is.null(control$maxfeval)) control$maxfeval<-2000*n
  if (is.null(control$eps)) control$epsl<-1e-07
  steptol <- control$eps
  # Hooke and Jeeves with bounds and masks
  if (length(upper) == 1) upper <- rep(upper, n)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (is.null(bdmsk)) { 
      idx <- 1:n 
  } else { idx <- which(bdmsk != 0) } # define masks
  nac <- length(idx)
  stepsize <- 1 # initial step size (could put in control())
  stepredn <- .1 # defined step reduction (again in control()??)
  offset = 100. # get from control() -- used for equality check
  if (any(par < lower) || any(par > upper)) stop("hjn: initial parameters out of bounds")
  pbase <- par # base parameter set (fold is function value)
  fold <- fn(par, ...) # fn at base point
  fcount <- 1 # count function evaluations, compare with maxfeval
  fmin <- fold # "best" function so far
  keepgoing <- TRUE
  ccode <- 1 # start assuming won't get to solution before feval limit
  while (keepgoing) {    
    # exploratory search
    for (jj in 1:nac){ # possibly could do this better in R
       # use unmasked parameters
       j <- idx[jj]
       vtemp <- par[j]
       donegstep <- FALSE
       if (vtemp + offset < upper[j] + offset) { # Not on upper bound 
          par[j] <- min(vtemp+stepsize, upper[j])
          if (par[j] + offset != vtemp + offset) {
             fcount <- fcount + 1
             f <- fn(par, ...)
          } else { f <- fold } # to force neg step
          if (f >= fold) { donegstep <- TRUE }
       }  else { donegstep <- TRUE } # on upper bound
       if (donegstep) {
          if (vtemp + offset > lower[j] + offset) { # can do negative step
             par[j] <- max(vtemp - stepsize, lower[j])
             if (par[j] + offset != vtemp + offset) {
               fcount <- fcount + 1
               f <- fn(par, ...)
             } else { f <- fold }
          }
       } # end donegstep
       if (f < fmin) { # improved function value
         fmin <- f
         changed <- TRUE
       } else { par[j] <- vtemp } # restore parameter value, no change
    } # end loop on axes
    if (control$trace > 0) { 
       cat("axial search with stepsize =",stepsize,"  fn value = ",fmin,"  after ",fcount,"\n")
    }
    if (fmin < fold) { # success -- do pattern move
       ptemp <- par
       par <- (par - pbase) + par # ignores masks
       for (j in 1:n) { # could do bounds check in vector?
          if (par[j] > upper[j]) par[j] <- upper[j]
          if (par[j] < lower[j]) par[j] <- lower[j]
       }
       pbase <- par
       par <- ptemp
       fold <- fmin
    } else { # no success in Axial Seart, so reduce stepsize
       if (fcount >= control$maxfeval) {
          keepgoing <- FALSE # too many function evaluations
          break
       }
       # first check if point changed
       samepoint <- any((par + offset) == (pbase + offset))
       if (samepoint) { 
          stepsize <- stepsize*stepredn
          if (stepsize <= steptol) keepgoing <- FALSE
          ccode <- 0 # successful convergence
       } else { # return to old base point
          par <- pbase
       }
    }
  } # end keepgoing loop
  ans <- list(par=par, value=fmin, counts=c(fcount, NA), convergence=ccode)
}
