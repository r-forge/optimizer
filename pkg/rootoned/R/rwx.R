rwx <- function(fn=fn, gr=NULL, ri=NULL, method="uniroot", ftrace=TRUE, ...){
  ## A wrapper for a variety of rootfinders
  ## We will try to document here
  # fn = the function for which root is to be found
  # ri = interval (2 element vector) of lower and upper bound to root
  #    OR  if second element is NA, then a guess to start the rootfinding
  # method = character name of the rootfinder (list below)
  # ftrace = logical TRUE if method is to be traced
  # ... = additional arguments to fn or FnTrace (e.g., ifn)
  ## ?? do we want to do a preliminary search? Esp. for GUESS methods
  
  ## OUTPUT
  #  root
  #  froot
  #  rtol
  #  niter (may be same as fncount)
  #  fncount
  #  method
  # ?? NOT INCLUDED grcount -- only newt1d uses gr now

  groot<-list(ifn=0, igr=0, ftrace=ftrace, fn=fn, gr=gr, label=method)
  envroot <- list2env(groot)

FnTrace <- function(x, ...) { 
  # Substitute function to call when rootfinding
  # Evaluate fn(x, ...)
    val <- envroot$fn(x, ...)
    envroot$ifn <- envroot$ifn + 1 # probably more efficient ways
    if (envroot$ftrace) {
       cat("f(",x,")=",val," after ",envroot$ifn," ",envroot$label,"\n")
    }
    val
}

grTrace <- function(x, ...) { 
  # Substitute function to call when rootfinding
  # Evaluate fn(x, ...)
  val <- envroot$gr(x, ...)
  envroot$igr <- envroot$igr + 1 # probably more efficient ways
  if (envroot$ftrace) {
    cat("gr(",x,")=",val," after ",envroot$igr," ",envroot$label,"\n")
  }
  val
}

  ## Preliminary checks
  if (length(ri) != 2) {stop("ri - the interval in which root to be found, must be 2 elements")}
  if (is.na(ri[2])) {
      fguess = ri[1]
  } else if (ri[2] <= ri[1]){stop("Lower bound must be > upper bound")}

  
## Beginning of methods
  
  if (method == "uniroot"){
    tst <- uniroot(FnTrace, interval=ri, ...) # may still need dotargs for tolerances etc.
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$init.it <- NULL
    #   return(tst)
  }
  
  if (method == "root1d"){
    tst <- root1d(f=FnTrace, interval=ri, trace=ftrace, ...)
    tst$iter <- tst$fcount
    tst$fcount <- NULL
  }
  
  if (method == "zeroin"){
    tst <- zeroin(f=FnTrace, ri, trace=ftrace, ...)
  }
  
  if (method == "newt1d"){
    if (is.null(gr)) stop("newt1d in rootwrap MUST have gradient")
    fguess<-ri[1]
    if (ftrace) cat("fguess = ", fguess,"\n")
    tst <- newt1d(fn=FnTrace, gr=grTrace, x0=fguess, trace=ftrace, ...)
    tst$iter <- tst$itn
    tst$itn <- NULL
    tst$rtol <- NA
    #   return(tst)
  }
  
  if (method == "newton"){
    fguess<-ri[1]
    tst <- newton(FnTrace, fguess, ...) # newton has dot args in pracma
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
  }
  
  if (method == "bisect"){
    tst <- bisect(f=FnTrace, ri[1], ri[2], ...) # No dotargs in pracma
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    #   return(tst)
  }
  
  if (method == "secant"){
    if (ftrace) print(ri)
    fguess <- ri[1]
    guess2 <- (fguess+0.01*(abs(fguess)+1))
    if (ftrace) cat("Start secant with fguess=", fguess," guess2=",guess2,"\n")
    tst <- secant(f=FnTrace, fguess, guess2, ...) 
    # has dotargs in pracma BUT DOES NOT USE THEM
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    #   return(tst)
  }
  
  if (method == "regulaFalsi"){
    tst <- regulaFalsi(f=FnTrace, ri[1], ri[2], ...) # no dotargs in pracma
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
    #   return(tst)
  } 

   if (method == "muller"){
    tst <- muller(f=FnTrace, ri[1], ri[2], 0.5*(ri[1]+ri[2]), ...) # no dotargs in pracma
    tst$iter <- tst$niter
    tst$niter <- NULL
    tst$froot <- tst$fval
    tst$fval <- NULL
    tst$rtol <- tst$reltol
    tst$reltol <- NULL
  #   return(tst)
  } 
  
  
  if (method == "brent"){
    tst <- brent(f=FnTrace, ri[1], ri[2], ...) # no dotargs in pracma
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$iter <- tst$f.calls
    tst$f.calls <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
  } 

  if (method == "ridders"){
    tst <- ridders(f=FnTrace, ri[1], ri[2], ...) # no dotargs in pracma
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
  } 
  
  
  
## End methods
  tst$fncount <- envroot$ifn # save value
  tst$method <- method # save method
#  rm(envroot) # clear the local environment
  tst

}

