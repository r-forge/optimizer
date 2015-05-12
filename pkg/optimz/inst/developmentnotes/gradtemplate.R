  # gr MUST be provided
  if (is.null(gr)) {  # if gr function is not provided STOP 
    stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminu") 
  }
  if (is.character(gr)) { # assume numerical gradient
  # Convert string to function call, assuming it is a numerical gradient function
    if (ctrl$trace > 0) cat("WARNING: using gradient approximation '",gr,"'\n")
    mygr<-function(par=par, userfn=fn, ...){
        do.call(gr, list(par, userfn, ...))
    }
  } else { 
    mygr<-gr 
  } # end else