
bench.GNE.nseq <- function(xinit, ..., echo=FALSE, control=list())
{
	
	methods <- c("Newton", "Broyden")
	globals <- c("none", "gline", "qline", "pwldog", "dbldog")
	globnames <- c("pure", "geom. line search", "quad. line search", "Powell trust region", "Dbl. trust region")
	
	fullnames1 <- paste(rep(methods, each=length(globals)), rep(globnames, length(methods)), sep=" - ")
	
	reslist <- list()
	
	for(m in methods)
	for(g in globals)
	{
		if(echo)
			cat("Begin", paste(m,g), "\n")
		mytime <- system.time(
							  res <- GNE.nseq(xinit, ..., method=m, global=g, control=control) 
							  , gcFirst = TRUE)
		res <- c(res, method=m, global=g, time= mytime[3], name=paste(m,g))
		if(echo)
			cat("End\n")
		reslist <- c(reslist, list(res))
	}	
	
	method <- "Levenberg-Marquardt"
	LM.Params <- c("merit", "jacmerit", "min")
	globals <- c("none", "gline", "qline")
	globnames <- c("pure", "geom. line search", "quad. line search")
	
	fullnames2 <- paste("LM", paste(rep(LM.Params, each=length(globals)), rep(globnames, length(LM.Params)), sep=" - "))
		
	for(l in LM.Params)
	for(g in globals)
	{
		if(echo)
			cat("Begin", paste(m,g), "\n")
		
		mytime <- system.time(
							  res <- GNE.nseq(xinit, ..., method=method, 
										  control=c(control, list(LM.param=l)), global=g) 
							  , gcFirst = TRUE)
		res <- c(res, method=method, global=g, LM.param=l, time= mytime[3], name=paste(g,l))
		if(echo)
			cat("End\n")
		reslist <- c(reslist, list(res))
	}	
	
	fctcall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[1] )
	jaccall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[2] )
	x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
	normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
	comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
	checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )	
	codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )
	localmethods <- sapply(1:length(reslist), function(i) reslist[[i]]$method )
	globalmethods <- sapply(1:length(reslist), function(i) reslist[[i]]$global )
	
#	list(compres = data.frame( method= c(fullnames1, fullnames2), fctcall, jaccall, 
#							  comptime, normFx, codes, localmethods, globalmethods, 
#							  x ), reslist = reslist)
	
	data.frame( method= c(fullnames1, fullnames2), fctcall, jaccall, comptime, normFx, codes, 
			   localmethods, globalmethods, x )
}



bench.GNE.ceq <- function(xinit, ..., echo=FALSE, control=list())
{
	
	methods <- "PR"
	globals <- "gline"
	globnames <- "geom. line search"
	
	fullnames1 <- paste(rep(methods, each=length(globals)), rep(globnames, length(methods)), sep=" - ")
	
	reslist <- list()
	
	for(m in methods)
	for(g in globals)
	{
		if(echo)
			cat("Begin", paste(m,g), "\n")
		mytime <- system.time(
							  res <- GNE.ceq(xinit, ..., method=m, global=g, control=control) 
							  , gcFirst = TRUE)
		res <- c(res, method=m, global=g, time= mytime[3], name=paste(m,g))
		if(echo)
			cat("End\n")
		reslist <- c(reslist, list(res))
	}	
	
	methods <- "AS"
	globals <- "pwldog"
	globnames <- "Powell trust region"
	
	fullnames2 <- paste(rep(methods, each=length(globals)), rep(globnames, length(methods)), sep=" - ")
	
	for(m in methods)
	for(g in globals)
	{
		if(echo)
			cat("Begin", paste(m,g), "\n")
		
		mytime <- system.time(
							  res <- GNE.ceq(xinit, ..., method=m, control=control, global=g) 
							  , gcFirst = TRUE)
		res <- c(res, method=m, global=g, time= mytime[3], name=paste(m,g))
		if(echo)
			cat("End\n")
		reslist <- c(reslist, list(res))
	}	
	
	fctcall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[1] )
	jaccall <- sapply(1:length(reslist), function(i) reslist[[i]]$counts[2] )
	x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
	normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
	comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
	checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )	
	codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )
	localmethods <- sapply(1:length(reslist), function(i) reslist[[i]]$method )
	globalmethods <- sapply(1:length(reslist), function(i) reslist[[i]]$global )
	
	
	data.frame( method= c(fullnames1, fullnames2), fctcall, jaccall, comptime, normFx, codes, 
			   localmethods, globalmethods, x )
}





