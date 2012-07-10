#algorithm for minimum finding
minGap <- function(xinit, Gap, gradGap, method=c("BB","BFGS"), 
	control=list(), ...)
{
	cat("Depreciated function.\n")
	
	method <- match.arg(method, c("BB","BFGS"))
	
	#default control parameters
	con <- list(tol=1e-6, maxit=100, echo=TRUE, stepinit=1)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	if(missing(Gap) || missing(xinit) || missing(gradGap))
		stop("missing parameters.")
		
	myGNE <- try( minfinding(xinit, Gap, gradGap, method, control=con, ...) )
	
	if(class(myGNE) == "try-error")
	{
		cat("\n")	
		stop("BB algorithm for minimisation of the (regularized) gap function does not converge.")
	}else
		myGNE
	
}



