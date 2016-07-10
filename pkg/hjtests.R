rm(list=ls())

# hjn.R -- R implementation of J Nash BASIC HJG.BAS 20160705
hjn <- function(par, fn, lower=-Inf, upper=Inf, bdmsk=NULL, control=list(trace=0), ...){
  n <- length(par) # number of parameters
  if (is.null(control$trace)) control$trace <- 0 # just in case
  if (is.null(control$stepsize)) {
     stepsize <- 1 # initial step size (could put in control())
  } else { stepsize <- control$stepsize }
  # Want stepsize positive or bounds get messed up
  if (is.null(control$stepredn)) {
     stepredn <- .1 # defined step reduction (again in control()??)
  } else { stepredn <- control$stepredn }
  if (is.null(control$maxfeval)) control$maxfeval<-2000*n
  if (is.null(control$eps)) control$epsl<-1e-07
  steptol <- control$eps
  # Hooke and Jeeves with bounds and masks
  if (length(upper) == 1) upper <- rep(upper, n)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (is.null(bdmsk)) { 
      bdmsk <- rep(1,n)
      idx <- 1:n 
  } else { idx <- which(bdmsk != 0) } # define masks
  if (any(lower >= upper)){
      warning("hjn: lower >= upper for some parameters -- set masks")
      bdmsk[which(lower >= upper)] <- 0
      idx <- which(bdmsk != 0)
  }
  cat("bdmsk:")
  print(bdmsk)
  cat("idx:")
  print(idx)
  nac <- length(idx)
  offset = 100. # get from control() -- used for equality check
  if (any(par < lower) || any(par > upper)) stop("hjn: initial parameters out of bounds")
  pbase <- par # base parameter set (fold is function value)
  f <- fn(par, ...) # fn at base point
  fmin <- fold <- f # "best" function so far
  pbest <- par # Not really needed 
  fcount <- 1 # count function evaluations, compare with maxfeval
    cat(fcount, "  f=",fold," at ")
    print(par)
#    tmp <- readline("cont.")
  keepgoing <- TRUE
  ccode <- 1 # start assuming won't get to solution before feval limit
  while (keepgoing) {    
    # exploratory search -- fold stays same for whole set of axes
    if (control$trace > 0) cat("Exploratory move - stepsize = ",stepsize,"\n")
    if (control$trace > 1) {
       cat("p-start:")
       print(par)
    }
    for (jj in 1:nac){ # possibly could do this better in R
       # use unmasked parameters
       j <- idx[jj]
       ptmp <- par[j]
       doneg <- TRUE # assume we will do negative step
       if (ptmp + offset < upper[j] + offset) { # Not on upper bound so do pos step 
          par[j] <- min(ptmp+stepsize, upper[j])
          if ((par[j] + offset) != (ptmp + offset)) {
             fcount <- fcount + 1
             f <- fn(par, ...)
               cat(fcount, "  f=",f," at ")
               print(par)
             if (f < fmin) {
                fmin <- f
                pbest <- par
                  cat("*")
                doneg <- FALSE # only case where we don't do neg
                resetpar <- FALSE
             } 
#             tmp <- readline("cont>")
          } 
       } # end not on upper bound
       if (doneg) {
         resetpar <- TRUE # going to reset the paramter unless we improve
         if ((ptmp + offset) > (lower[j] + offset)) { # can do negative step
            par[j] <- max((ptmp - stepsize), lower[j])
            if ((par[j] + offset) != (ptmp + offset)) {
               fcount <- fcount + 1
               f <- fn(par, ...)
                 cat(fcount, "  f=",f," at ")
                 print(par)
               if (f < fmin) {
                  fmin <- f
                  pbest <- par
                  cat("*")
                  resetpar <- FALSE # don't reset parameter
               } 
#              tmp <- readline("cont<")
            }
         } #  neg step possible
       } # end doneg
       if (resetpar) { par[j] <- ptmp }
    } # end loop on axes
    if (control$trace > 0) { 
       cat("axial search with stepsize =",stepsize,"  fn value = ",fmin,"  after ",fcount,"\n")
    }
    if (fmin < fold) { # success -- do pattern move (control$trace > 0) cat("Pattern move \n")
       if (control$trace > 1) {
          cat("PM from:")
          print(par)
          cat("pbest:")
          print(pbest)
       }
       for (jj in 1:nac) { # Note par is best set of parameters
          j <- idx[jj]
          ptmp <- 2.0*par[j] - pbase[j]
          if (ptmp > upper[j]) ptmp <- upper[j]
          if (ptmp < lower[j]) ptmp <- lower[j]
          pbase[j] <- par[j]
          par[j] <- ptmp 
       }
       fold <- fmin
       if (control$trace > 1) {
          cat("PM to:")
          print(par)
       }
    # Addition to HJ -- test new  base
#       fcount <- fcount + 1
#       f <- fn(par, ...)
#         cat(fcount, "  f=",f," at ")
#         print(par)
#         tmp <- readline("PM point")
#       if (f < fmin) {
#         if (control$trace > 0) {cat("Use PM point as new base\n")}
#         pbest <- pbase <- par
#       }
    } else { # no success in Axial Seart, so reduce stepsize
       if (fcount >= control$maxfeval) {
          keepgoing <- FALSE # too many function evaluations
          break
       }
       # first check if point changed
       samepoint <- identical((par + offset),(pbase + offset))
       if (samepoint) { 
          stepsize <- stepsize*stepredn
          cat("Reducing step to ",stepsize,"\n")
          if (stepsize <= steptol) keepgoing <- FALSE
          ccode <- 0 # successful convergence
       } else { # return to old base point
          if (control$trace > 1) {
             cat("return to base at:")
             print(pbase)
          }
          par <- pbase
       }
    }
  } # end keepgoing loop 
  if (identical(pbest, pbase)) {cat("pbase = pbest\n") }
  else { cat("BAD!: pbase != pbest\n") }
   
  ans <- list(par=pbest, value=fmin, counts=c(fcount, NA), convergence=ccode)
}

fn <- function(x) sum(x*x)

ans1 <- hjn(c(1, 1), fn, control=list(trace=1))

print(ans1)

library(adagio)   # fnRosenbrock

ans2 <- hjn(rep(0,10), fnRosenbrock, lower = rep(-5.12,10), upper = rep(5.12,10), control=list(trace=2, stepredn=0.2))
print(ans2)

library(dfoptim)
ans2h <- hjkb(rep(0,10), fnRosenbrock, lower = rep(-5.12,10), upper = rep(5.12,10), control=list(info=2))
print(ans2h)

library(dfoptim)  # hjkb

ans3 <- hjkb(rep(0, 10), fnRosenbrock, lower = rep(-5.12,10), upper = rep(5.12,10), control=list(info=1))

# Simple bounds test for n=4
bt.f<-function(x){
 sum(x*x)
}
bt.g<-function(x){
  gg<-2.0*x
}
n<-4
lower<-rep(0,n)
upper<-lower # to get arrays set
bdmsk<-rep(1,n)
# bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
    lower[i]<-1.0*(i-1)*(n-1)/n
    upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)

alhn<-hjn(xx, bt.f, lower=lower, upper=upper, control=list(trace=2))
alhn

cat("Now force a mask upper=lower for parameter 3 and see what happens\n")
upper[3] <- lower[3]
xx[3] <- lower[3] # and set parameter

am <- hjn(xx, bt.f, lower=lower, upper=upper, control=list(trace=2))
am

require(optimr)
allbdm <- opm(xx, bt.f, bt.g, lower=lower, upper=upper, method="ALL", control=list(trace=0))
summary(allbdm, order=value)



