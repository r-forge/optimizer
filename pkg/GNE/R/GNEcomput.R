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

GNE.nseq <- function(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), method="Newton", control=list(), ...)
{
	if(method == "default")
		method <- "Newton"
		
	nseq(xinit, Phi, jacPhi, argPhi, argjac, method, control, ...)	
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



bench.GNE.nseq <- function(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), echo=FALSE, ...)
{
	
methods <- c("Newton", "Broyden")
globals <- c("none", "gline", "qline", "pwldog", "dbldog")
globnames <- c("pure", "geom. line search", "quad. line search", "Powell trust region", "Dbl. trust region")

fullnames <- paste(rep(methods, each=length(globals)), rep(globnames, length(methods)), sep=" - ")


times <- vector("numeric", length(methods) * length(globals) )
reslist <- list()

for(m in methods)
for(g in globals)
{
mytime <- system.time(
	 res <- GNE.nseq(xinit, Phi, jacPhi, argPhi, argjac, ..., method=m, global=g) 
, gcFirst = TRUE)
	res <- c(res, method=m, global=g, time= mytime[3])
	if(echo)
		print(res)
	reslist <- c(reslist, list(res))
}	

fctcall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[1] )
jaccall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[2] )
x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )

list(compres = data.frame( method= fullnames, fctcall, jaccall, comptime, x,normFx ),
	reslist = reslist )
}

