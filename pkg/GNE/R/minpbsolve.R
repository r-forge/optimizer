#algorithm for minimum finding
minpb <- function(xinit, Gap, gradGap, argGap, argGrad, method=c("BB","BFGS"), 
	control=list(), ...)
{
	method <- match.arg(method, c("BB","BFGS"))
	
	#default control parameters
	con <- list(tol=1e-6, maxit=100, echo=TRUE, stepinit=1)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	if(missing(Gap) || missing(xinit) || missing(gradGap))
		stop("missing parameters.")

	
	wrapGap <- function(x, argGap, argGrad)
		evalwitharglist(Gap, x, argGap) 
	
	wrapGrad <- function(x, argGap, argGrad)
		evalwitharglist(gradGap, x, argGrad)
	
		
	myGNE <- try( minfinding(xinit, wrapGap, wrapGrad, method, control=con, ...),
				silent = TRUE)
	
	if(class(myGNE) == "try-error")
	{
		myGNE <- "BB algorithm for minimisation of the (regularized) gap function does not converge."
	}else
	{
		myGNE
	}
}



