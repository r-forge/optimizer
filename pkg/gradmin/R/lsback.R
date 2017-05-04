lsback<-function(w, ...) { 
  # lsback -- backtrack line search
  st <- 1.0
  w$gproj <- as.numeric(crossprod(w$grd, w$tdir) ) # can we save this from elsewhere?
  if(w$trace > 1) cat("lsback: Gradient projection = ",w$gproj,"\n")
  cat("offset:",w$offset," tdir, xb\n")
  print(w$tdir)
  print(w$xb)
  repeat {
    w$xnew <- w$xb + st*w$tdir # new point
    if (all((w$offset+w$xnew) == (w$offset+w$xb))) { # no better parameters
      if (w$trace > 1) cat("No progress in lsback\n")
      st <- 0
      rlout <- st
      attr(rlout,"Fval")<-w$fbest # Assume we pass this in
      attr(rlout,"FAIL")<- TRUE
      return(rlout)
    }
    fval <- w$fn(w$xnew, ...)
    w$nf <- w$nf + 1
    if (w$trace > 2) cat("Step = ",st," fval = ",fval,"\n")
    if ((w$trace > 1) && (w$watch)) tmp <- readline("lnsrch continue?")
    if (fval <= w$fbest + w$acctol*st*w$gproj) break # Armijo condition
    st <- w$stepdec*st # new step
  }
  w$stp <- st # save this for working (??simplify!)
  rlout <- st
  attr(rlout, "Fval")<- fval # could return through w
  attr(rlout,"FAIL")<- FALSE
  rlout
} # end backtrack line search
