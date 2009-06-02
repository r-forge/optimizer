BBsolve <- function(par, fn, algorithm = c("dfsane", "sane"), 
        method=c(2,3,1), M = c(10, 50), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500) 
    {
    control.pars <- expand.grid(method=method, M=M, NM=NM)
    alg <- match.arg(algorithm)

    feval <- iter <-  0
    ans.best.value <- Inf
    for (i in 1: nrow(control.pars) ) {
      cpars <- unlist(control.pars[i, ])
      #cat("Try : ", i, "Method = ", cpars[1], "M = ", cpars[2], "Nelder-Mead = ", cpars[3], "\n")

      if (alg == "dfsane") temp <- 
 	dfsane(par=par, fn, method=cpars[1], control=list(M=as.numeric(cpars[2]), 
 	   NM=cpars[3], 
 	   maxit=maxit, tol=tol, trace=FALSE, noimp=min(100, 5*cpars[2])), ...)

      if (alg == "sane") temp <-
 	  sane(par=par, fn, method=cpars[1], control=list(M=as.numeric(cpars[2]),
 	   NM=cpars[3], 
 	   maxit=maxit, tol=tol, trace=FALSE, noimp=min(100, 5*cpars[2])), ...)

      feval <- feval + temp$feval
      iter <- iter + temp$iter

      if (temp$convergence  == 0) {
 	   ans.best <- temp
 	   ans.best$feval <- feval
 	   ans.best$iter <- iter
 	   ans.best$cpar <- cpars
 	   break
 	   } 
      else if (temp$residual < ans.best.value) {
 	   ans.best <- temp
 	   ans.best.value <- ans.best$residual
 	   ans.best$feval <- feval
 	   ans.best$iter <- iter
 	   ans.best$cpar <- cpars
 	   }

      }  # "i" loop completed

    if (ans.best$convergence != 0)
         cat ("  Unsuccessful convergence.\n")
    else cat ("  Successful convergence.\n")

    ans.best
    }

