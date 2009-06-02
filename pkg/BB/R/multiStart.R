multiStart <- function(par, fn, algorithm = c("dfsane", "sane"), 
	method=c(2,3,1), M = c(10, 50), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500, details=FALSE) 
    {
     if (is.null(dim(par))) par <- matrix(par, nrow=1, ncol=length(par))
    if (ncol(par) > 20) NM <- FALSE   
    
    ans <- vector("list",    length=nrow(par))
    cvg <- vector("logical", length=nrow(par))
    values <- vector("numeric", length=nrow(par))
    
    feval <- iter <-  0
    pmat <- matrix(NA, nrow(par) ,ncol(par))
    
    for (k in 1:nrow(par)){
       cat("Parameter set : ", k, "... \n")
    
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
    
    #   some duplication here (everything is in info=ans)
    if (details) list(par=pmat, values=values, converged=cvg, info=ans)
    else list(par=pmat, values=values, converged=cvg)
    }
