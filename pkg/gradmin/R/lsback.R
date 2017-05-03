lsback<-function(fn, fbest, xc, d, grv, w, ...) { 
  # lsback -- backtrack line search
  st <- 1.0
  gproj <- as.numeric(crossprod(grv, xc) )
  repeat {
    xnew <- xc + st*d # new point
    if (all((w$offset+xnew) == (w$offset+xc))) { # no better parameters
      st <- 0
      rlout <- st
      attr(rlout,"Fval")<-fbest # Assume we pass this in
      attr(rlout,"FAIL")<- TRUE
      return(rlout)
    }
    fval <- fn(xnew, ...)
    w$nf <- w$nf + 1
    if (w$trace > 1) cat("Step = ",st," fval = ",fval,"\n")
    if (fval <= fbest + w$acctol*st*gproj) break # Armijo condition
    st <- w$stepdec*st # new step
  }
  rlout <- st
  attr(rlout, "Fval")<- fval
  attr(rlout,"FAIL")<- FALSE
  rlout
} # end backtrack line search
