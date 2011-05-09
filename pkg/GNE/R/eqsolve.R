eqsolve <- function(xinit, f, jac, control=list())
{
	#default control parameters
	con <- list(tol=1e-6, maxit=100, echo=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	xk <- xinit
	fk <- f(xinit)
	Jacfk <- jac(xinit)
	iter <- 0
	
	
	while(sqrt( sum( fk^2 ) ) > con$tol && iter < con$maxit) 
	{
		xk <- NewtonNext(xk, fk, Jacfk)
	
		fk <- f(xk)
		Jacfk <- jac(xk)
		
		iter <- iter+1
		if(con$echo >= 1)
			cat("**** k", iter, "\n x_k", xk, "\n")
		if(con$echo >= 2)
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
		if(con$echo >= 3)
		{
			cat(" Jac Phi(x_k)\n")
			print(Jacfk)
		}
	

	}

	list(par= as.vector(xk), value=sqrt( sum( fk^2 ) ), counts=iter, iter=iter, 
		 code=(iter >= con$maxit)*1, message=NULL)

}