gridStart <- function(par, fn, algorithm = c("dfsane", "sane"), 
	method=c(2,1,3), M = c(10, 100), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500) 
    {
    par <- as.matrix(par)
    if (ncol(par) > 20) NM <- FALSE   
    
    ans <- vector("list",    length=nrow(par))
    cvg <- vector("logical", length=nrow(par))
    values <- vector("numeric", length=nrow(par))
    
    feval <- iter <-  0
    pmat <- matrix(NA, nrow(par) ,ncol(par))
    
    for (k in 1:nrow(par)){
       cat("Parameter set : ", k, "... ")
    
       ans[[k]] <- BBsolve(par[k,], fn, algorithm = algorithm, 
 	    method=method, M = M, NM=NM, ... , tol=tol, maxit=maxit) 
    
       cvg[k]  <-  (ans[[k]]$convergence  == 0)
       values[k] <- ans[[k]]$residual
       pmat[k, ] <- ans[[k]]$par
       }  # "k" loop completed

    #  need to think about what should really be returned here
    #   possibly some stats about the distribution of solutions
    #     around ans[[which.min(values)]]
    #  if (!any(cvg)) ans[[1]] <- ans.best
    
    #success <- !is.na(pmat[,1])  # this was indicating convergence (cvg)
      
    #  #  I don't understand k==1 here (from old code).  I think
    #  #  this condition may only happen when there was a single par setting
    #  #  which would now be done with a direct call to BBsolve, not gridStart.
    #  if (k ==1 | !any(cvg)) ans[[1]] 
    #  else list(par = pmat[success, ], info=ans[success])
    
    # # if (!any(cvg))  ans.best else
    #list(par=pmat[cvg,], info=ans[cvg]) #think I want failures too

    #   some duplication here (everything is in info=ans)
    list(par=pmat, values=values, converged=cvg, info=ans)
    }



BBsolve <- function(par, fn, algorithm = c("dfsane", "sane"), 
	method=c(2,1,3), M = c(10, 100), NM=c(TRUE, FALSE), ... , 
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

