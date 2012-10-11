#algorithm for minimum finding
minpb <- function(xinit, fn, argfn, gr=NULL, arggr, 
	hin=NULL, arghin, hin.jac=NULL, arghin.jac, method=c("BB","BFGS"), 
	control=list(), silent=TRUE, ...)
{
	method <- match.arg(method, c("BB","BFGS"))
	if(missing(argfn))
		argfn <- NULL
	if(missing(arggr))
		arggr <- NULL
	if(missing(arghin))
		arghin <- NULL
	if(missing(arghin.jac))
		arghin.jac <- NULL
	
	#default control parameters
	conBB <- list(tol=1e-8, maxit=100, echo=FALSE, stepinit=1)
	namcBB <- names(conBB)
	conBB[namcBB <- names(control)] <- control
	conoptim <- list(abstol=1e-8, maxit=100, trace=0)
	namcoptim <- names(conoptim)
	conoptim[namcoptim <- names(control)] <- control
	conalabama <- list(eps=1e-8, itmax=100, trace=FALSE)
	namcalabama <- names(conalabama)
	conalabama[namcalabama <- names(control)] <- control

	if(missing(fn) || missing(xinit))
		stop("missing parameters.")

	if(method == "BB")
		stop("not yet implemented.")

	#create the objective function
	Fn <- function(x) evalwitharglist(fn, x, argfn)
	test.try <- try( Fn(xinit), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate fn(init).", fvec=NA) )
	if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="fn(init) has infinite, NA or NaN values.", fvec=NA) )
	if(!silent)
		cat("test of function done.\n")
	
	#create the gradient function (if provided)
	if(!is.null(gr))
		Gr <- function(x) evalwitharglist(gr, x, arggr)
	else
		Gr <- NULL
	
#	print(attributes(arghin))
#	print(evalwitharglist(hin, xinit, arghin))
	#create the inequality constraint
	if(!is.null(hin))
	{
		Hin <- function(x) evalwitharglist(hin, x, arghin)
		test.try <- try( Hin(xinit), silent=silent )
		if(class(test.try) == "try-error")
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="Can't evalate hin(init).", fvec=NA) )
		if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="hin(init) has infinite, NA or NaN values.", fvec=NA) )
	}else
		Hin <- NULL
	if(!silent)
		cat("test of inequality constraint done.\n")
	
	#create the inequality constraint Jacobian
	if(!is.null(hin.jac))
		Hin.jac <- function(x) evalwitharglist(hin.jac, x, arghin.jac)
	else
		Hin.jac <- NULL

	if(!silent)
		cat("beginning of optimization:\t")
	
	if(is.null(Hin))
	{
		if(!silent)
			cat("optimization carried out by optim\n")
		test.try <- try( optim(par=xinit, fn=Fn, gr=Gr, control=conoptim, 
							   method=method), silent=silent)
		if(class(test.try) == "try-error")
			res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."))
		if(class(test.try) != "try-error")
			res <- list(par = test.try$par, value = test.try$value, 
					counts = test.try$counts, iter = as.numeric(test.try$counts["function"]), 
					code = test.try$convergence, message = test.try$message)		
	}else 
	{
		require(alabama)
		if(!silent)
			cat("optimization carried out by constrOptim.nl\n")
		test.try <- try( constrOptim.nl(par=xinit, fn=Fn, gr=Gr, hin=Hin, 
						hin.jac=Hin.jac, control.outer=conalabama), silent=silent)
		if(class(test.try) == "try-error")
			res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."))
		if(class(test.try) != "try-error")
			res <- list(par = test.try$par, value = test.try$value, 
					counts = test.try$counts, iter = test.try$outer.iterations, 
					code = test.try$convergence, message = test.try$message)		
		
	}
	if(!silent)
		cat("end of optimization\n")
	
	res
}



