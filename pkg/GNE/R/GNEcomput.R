GNE <- function(approach = c("non smooth", "fixed point", "minimization"), method = "default", xinit,
	func1, func2, argfunc1=list(), argfunc2=list(), ..., control=list())
{
	approach <- match.arg(approach, c("non smooth", "fixed point", "minimization"))
	
	if(approach == "non smooth")
		res <- GNE.nseq(xinit, func1, func2, argfunc1, argfunc2, method, control, ...)
	
	if(approach == "fixed point")
		res <- GNE.fpeq(xinit, func1, func2, argfunc1, argfunc2, method, control, ...)
	
	if(approach == "minimization")
		res <- GNE.min(xinit, func1, func2, argfunc1, argfunc2, method, control, ...)
	
	c(res, approach=approach)
}

GNE.nseq <- function(xinit, phi, jacphi, argphi=list(), argjac=list(), method="Newton", control=list(), ...)
{
	if(method == "default")
		method <- "Newton"
		
	nseq(xinit, phi, jacphi, argphi, argjac, method, control, ...)	
}


GNE.fpeq <- function(xinit, yfunc, func2, argy=list(), argfunc2=list(), method, control=list(), ...)
{
	stop("not yet implemented.")
	if(method == "default")
		method <- "MPE"
		
	fpeq(xinit, yfunc, func2, argy, argfunc2, method, control, ...)
}


GNE.min <- function(xinit, gap, gradgap, arggap=list(), arggrad=list(), method, control=list(), ...)
{
	if(method == "default")
		method <- "BB"
	
	minpb(xinit, gap, gradgap, arggap, arggrad, method, control, ...)
}
