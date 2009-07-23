BBoptim <- function(par, fn, gr=NULL, method=c(2,3,1), project=NULL, 
     lower=-Inf, upper=Inf, control=list(), quiet=FALSE, ...) 
    {
    ctrl <- list(maxit = 1500, M = c(50, 10), gtol = 1e-05, maxfeval = 10000, 
        maximize = FALSE, trace = FALSE, triter = 10, eps = 1e-07)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if(is.matrix(par)) stop("argument par should not be a matrix in BBoptim.")
    ctrl[namc] <- control
    M <- ctrl$M
    maxit <- ctrl$maxit
    gtol <- ctrl$gtol
    maxfeval <- ctrl$maxfeval
    maximize <- ctrl$maximize
    trace <- ctrl$trace
    triter <- ctrl$triter
    eps <- ctrl$eps
    control.pars <- expand.grid(method=method, M=M)

    feval <- iter <-  0
    ans.best.value <- Inf
    for (i in 1: nrow(control.pars) ) {
      cpars <- unlist(control.pars[i, ])
      temp <- try(spg(par=par, fn=fn, gr=gr, method=cpars[1], project=project, 
	            lower=lower, upper=upper, 
		    control=list(M=as.numeric(cpars[2]), maxit=maxit, 
		       maximize=maximize, trace=trace, triter=triter, 
		       maxfeval=maxfeval, eps=eps), ...),    silent=TRUE)

      if (!inherits(temp, "try-error")) {
   	 feval <- feval + temp$feval
   	 iter <- iter + temp$iter

   	 if (temp$convergence  == 0) {
   	      ans.best <- temp
   	      ans.best$feval <- feval
   	      ans.best$iter <- iter
   	      ans.best$cpar <- cpars
   	      break
   	      } 
   	 else if (temp$value < ans.best.value) {
   	      ans.best <- temp
   	      ans.best.value <- ans.best$value
   	      ans.best$feval <- feval
   	      ans.best$iter <- iter
   	      ans.best$cpar <- cpars
   	      }
         }
      }  # "i" loop completed

    if(!quiet) {if (ans.best$convergence != 0)
                     cat ("  Unsuccessful convergence.\n")
                else cat ("  Successful convergence.\n")
		}

    ans.best
    }