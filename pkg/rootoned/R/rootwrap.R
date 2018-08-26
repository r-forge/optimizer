rootwrap <- function(fn=fn, gr=NULL, ri=NULL, method="uniroot", ftrace=TRUE, ...){
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
  #  grcount
  
  ## Preliminary checks
  # print(ri)
  if (length(ri) != 2) {stop("ri - the interval in which root to be found, must be 2 elements")}
  if (is.na(ri[2])) {
      fguess = ri[1]
  } else if (ri[2] <= ri[1]){stop("Lower bound must be > upper bound")}
  TraceSetup(ftrace=ftrace) # turn on the trace
  
## Beginning of methods
  
  if (method == "uniroot"){
    envroot$fn <- fn
    envroot$label <- "uniroot"
    tst <- uniroot(FnTrace,ri,  ...)
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$init.it <- NULL
    #   return(tst)
  }
  
  if (method == "root1d"){
    if( ! require(rootoned, quietly=TRUE) ) stop("Package rootoned for root1d not found")
    envroot$fn <- fn
    envroot$label <- "root1d" 
    tst <- root1d(FnTrace,ri, trace=ftrace, ...)
    tst$iter <- tst$fcount
    tst$fcount <- NULL
    ## ?? be nice to include trace info from method like bisect and fp in FnTrace, but...
    #   return(tst)
  }
  
  if (method == "zeroin"){
    if( ! require(rootoned, quietly=TRUE) ) stop("Package rootoned for zeroin not found")
    envroot$fn <- fn
    envroot$label <- "rootoned::newt1d"
    tst <- zeroin(FnTrace, ri, trace=ftrace, ...)
    tst$iter <- tst$maxit
    tst$maxit <- NULL
    #   return(tst)
  }
  
  if (method == "newt1d"){
    if( ! require(rootoned, quietly=TRUE) ) stop("Package rootoned for root1d not found")
    envroot$fn <- fn
    envroot$label <- "rootoned::newt1d"
    # gr <- function(x) { cos(x)}
    fguess<-ri[1]
    tst <- newt1d(FnTrace, fguess, gr=gr, label="newt1d", trace=ftrace, ...)
    tst$iter <- tst$itn
    tst$itn <- NULL
    tst$rtol <- NA
    #   return(tst)
  }
  
  if (method == "newton"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for newton not found")
    envroot$fn <- fn
    envroot$label <- "pracma::newton"
    # gr <- function(x) { cos(x)} ?? Why does newton not need gr?
    fguess<-ri[1]
    tst <- newton(FnTrace, fguess, ...)
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
  }
  
  if (method == "bisect"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for bisect not found")
    envroot$fn <- fn
    envroot$label <- "pracma::bisect"
    tst <- bisect(f=FnTrace, ri[1], ri[2])
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    #   return(tst)
  }
  
  if (method == "secant"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for secant not found")
    envroot$fn <- fn
    envroot$label <- "pracma::secant"
    fguess <- ri[1]
    tst <- secant(f=FnTrace, fguess, (fguess+0.01*(abs(fguess)+1)), ...)
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    #   return(tst)
  }
  
  if (method == "regulaFalsi"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for regulaFalsi not found")
    envroot$fn <- fn
    envroot$label <- "pracma::regulaFalsi"
    tst <- regulaFalsi(f=FnTrace, ri[1], ri[2])
    tst$froot <- tst$f.root
    tst$f.root <- NULL
    tst$rtol <- tst$estim.prec
    tst$estim.prec <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
    #   return(tst)
  } 

   if (method == "muller"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for muller not found")
    envroot$fn <- fn
    envroot$label <- "pracma::muller"
    tst <- muller(f=FnTrace, ri[1], ri[2])
    tst$iter <- tst$niter
    tst$niter <- NULL
    tst$froot <- tst$fval
    tst$fval <- NULL
    tst$rtol <- tst$reltol
    tst$reltol <- NULL
  #   return(tst)
  } 
  
  
  if (method == "brent"){
    if( ! require(pracma, quietly=TRUE) ) stop("Package pracma for brent(Dekker) not found")
    envroot$fn <- fn
    envroot$label <- "pracma::brent"
    tst <- muller(f=FnTrace, ri[1], ri[2])
    tst$froot <- tst$fval
    tst$fval <- NULL
#    tst$froot <- tst$f.root
#    tst$f.root <- NULL
    tst$iter <- tst$f.calls
    tst$f.calls <- NULL
    tst$rtol <- tst$reltol
    tst$reltol <- NULL
#    tst$rtol <- tst$estim.prec
#    tst$estim.prec <- NULL
    tst$iter <- tst$niter
    tst$niter <- NULL
  #   return(tst)
  } 
  
  
## End methods
  tst$fncount <- envroot$ifn # save value
  tst$method <- method # save method
#  rm(envroot) # clear the local environment
  tst

}

