lsback<-function(w, ...) { 
  # lsback -- backtrack line search
  st <- 1.0
  w$gproj <- as.numeric(crossprod(w$grd, w$tdir) ) # can we save this from elsewhere?
#  if(w$trace > 1) cat("Gradient projection = ",w$gproj,"\n")
  repeat {
    xnew <- w$xb + st*w$tdir # new point
    if (all((w$offset+xnew) == (w$offset+w$xc))) { # no better parameters
      if (w$trace > 1) cat("No progress in lsback\n")
      st <- 0
      rlout <- st
      attr(rlout,"Fval")<-w$fbest # Assume we pass this in
      attr(rlout,"FAIL")<- TRUE
      return(rlout)
    }
    fval <- fn(xnew, ...)
    w$nf <- w$nf + 1
    if (w$trace > 2) cat("Step = ",st," fval = ",fval,"\n")
    if ((w$trace > 1) && (w$watch)) tmp <- readline("lnsrch continue?")
    if (fval <= w$fbest + w$acctol*st*gproj) break # Armijo condition
    st <- w$stepdec*st # new step
  }
  rlout <- st
  attr(rlout, "Fval")<- fval # could return through w
  attr(rlout,"FAIL")<- FALSE
  rlout
} # end backtrack line search
