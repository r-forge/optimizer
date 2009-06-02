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
