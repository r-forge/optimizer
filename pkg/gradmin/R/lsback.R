lsback<-function(fn, fbest, xc, d, grv, ws, ...) { 
  # lsback -- backtrack line search
  st <- 1.0
  gproj <- as.numeric(crossprod(grv, xc) )
  repeat {
    xnew <- xc + st*d # new point
    if ((ws$offset+xnew) == (ws$offset+xc)) { # no better parameters
      st <- 0
      rlout <- st
      attr(rlout,"Fval")<-fbest # Assume we pass this in
      return(rlout)
    }
    fval <- fn(xnew, ...)
    ws$nf <- ws$nf + 1
    if (ws$trace > 1) cat("Step = ",st," fval = ",fval,"\n")
    if (fval <= fbest + ws$acctol*st*gproj) break # Armijo condition
    st <- ws$stepdec*st # new step
  }
  rlout <- st
  attr(rlout, "Fval")<- fval
  rlout
} # end backtrack line search
