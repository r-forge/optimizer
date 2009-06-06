multiStart <- function(par, fn, gr=NULL, action = c("solve", "optimize"),
	method=c(2,3,1), control=list(), ... ) 
    {
     if (is.null(dim(par))) par <- matrix(par, nrow=1, ncol=length(par))
    
    ans <- vector("list",    length=nrow(par))
    cvg <- vector("logical", length=nrow(par))
    values <- vector("numeric", length=nrow(par))
    
    action <- match.arg(action)

    feval <- iter <-  0
    pmat <- matrix(NA, nrow(par) ,ncol(par))
    
    for (k in 1:nrow(par)){
       cat("Parameter set : ", k, "... \n")
    
       if (action == "solve") ans[[k]] <- BBsolve(par[k,], fn, method=method, 
		control=control, ...) 
       if (action == "optimize") ans[[k]] <- BBoptim(par[k,], fn=fn, gr=gr, method=method, 
		control=control, ...) 
    
       cvg[k]  <-  (ans[[k]]$convergence  == 0)
        values[k] <- if (action == "solve") ans[[k]]$residual
		                      else  ans[[k]]$value
       pmat[k, ] <- ans[[k]]$par
       }  # "k" loop completed

    ans.ret <- list(par=pmat, fvalue=values, converged=cvg)
    attr(ans.ret, "details") <- ans
    ans.ret
    }
