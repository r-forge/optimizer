GNE <- function(approach = 
	c("non smooth", "fixed point", "minimization", "constr equation"), 
	method = "default", xinit, func1, func2, argfunc1=list(), 
	argfunc2=list(), ..., control=list())
{
	approach <- match.arg(approach, c("non smooth", "fixed point", "minimization", "constr equation"))
	
	if(approach == "non smooth")
		res <- GNE.nseq(xinit, func1, func2, argfunc1, argfunc2, method, control, ...)

	if(approach == "constr equation")
		res <- GNE.ceq(xinit, argfunc1, argfunc2, method, control, ...)
	
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

GNE.ceq <- function(xinit, argH, argjacH, method="IP", control=list(), ...)
{
	if(method == "default")
		method <- "IP"
		
	ceq(xinit, argH, argjacH, method, control, global="gline", ...)	
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
	res <- c(res, method=m, global=g, time= mytime[3], name=paste(m,g))
	if(echo)
		print(res)
	reslist <- c(reslist, list(res))
}	

fctcall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[1] )
jaccall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[2] )
x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )	
codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )		

list(compres = data.frame( method= fullnames, fctcall, jaccall, comptime, x, normFx, codes ),
	reslist = reslist)
}


bench.GNE.nseq.LM <- function(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), 
	echo=FALSE, control=list())
{
	
	method <- "Levenberg-Marquardt"
	LM.Params <- "min"  #c("merit", "jacmerit", "min")
	globals <- c("none", "gline", "qline")
	globnames <- c("pure", "geom. line search", "quad. line search")
	
	fullnames <- paste("LM", paste(rep(LM.Params, each=length(globals)), rep(globnames, length(LM.Params)), sep=" - "))
		
	
	times <- vector("numeric", length(LM.Params) * length(globals) )
	reslist <- list()
	
	
	for(l in LM.Params)
	for(g in globals)
	{
		
		mytime <- system.time(
			res <- nseq(xinit, Phi, jacPhi, argPhi, argjac, method=method, 
						control=c(control, list(LM.param=l)), global=g) 
							  , gcFirst = TRUE)
		res <- c(res, method=method, global=g, LM.param=l, time= mytime[3], name=paste(g,l))
		if(echo)
			print(res)
		reslist <- c(reslist, list(res))
	}	
	mytime <- system.time(
		res <- nseq(xinit, Phi, jacPhi, argPhi, argjac, method=method, 
					  control=c(control, list(LM.param="adaptive")), global="none") 
						  , gcFirst = TRUE)
	res <- c(res, method=method, global="none", LM.param="adaptive", time= mytime[3], name=paste("none", "adaptive"))
	if(echo)
		print(res)
	reslist <- c(reslist, list(res))
	fullnames <- c(fullnames, paste("LM", paste("adaptive", "pure", sep=" - ")) )
	
	
	fctcall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[1] )
	jaccall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[2] )
	x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
	normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
	comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
	checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )
	codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )	
	
	list(compres = data.frame( method= fullnames, fctcall, jaccall, comptime, x, normFx, codes ),
		 reslist = reslist )
}

