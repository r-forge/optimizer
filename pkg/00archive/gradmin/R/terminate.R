terminate <- function(w, ...){
  if (w$trace > 2) cat("Termination test\n")    
  halt <- FALSE # default is keep going
  # tests on too many counts??
  if (w$niter >= w$maxit) {
    if (w$trace > 0) cat("Too many (",niter," iterations\n")
    halt <- TRUE
    w$convcode <- 1
    if (w$trace > 2) cat("returning X niter\n")
    return(halt)
  }
  if (w$nf >= w$maxfevals) {
    w$msg <- paste("Too many (",w$nf,") function evaluations")
    if (w$trace > 0) cat(w$msg,"\n")
    halt <- TRUE
    w$convcode <- 1 # ?? value
    if (w$trace > 2) cat("returning X nf\n")
    return(halt)
  }
  #    if (w$ng > w$maxgevals){} # not implemented
  #    if (w$nh > w$maxhevals){} # not implemented
  gmax <- max(abs(w$grd))
  if (gmax <= w$epstol*w$f0) {
    w$msg <- paste("Small gradient norm",gmax)
    if (w$trace > 0) cat(w$msg,"\n")
    halt <- TRUE
    w$convcode <- 0 # OK
    if (w$trace > 2) cat("returning X gmax\n")
    return(halt)
  }
  if (w$lsfail) {
    if (w$ng > w$lastsd) { 
      if (w$trace > 1) cat("lsfail but need to try steepest descents")
      w$resetB <- TRUE # go round again on steepest descent
    } else { 
      halt <-TRUE 
      if (w$trace > 1) cat("returning lsfail with ng, lastsd ",w$ng, w$lastsd,"\n")
    }
  }
  halt
}  
