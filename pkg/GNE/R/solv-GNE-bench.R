
bench.GNE.nseq <- function(xinit, ..., echo=FALSE, control=list())
{
	
	methods <- c("Newton", "Broyden")
	globals <- c("none", "gline", "qline", "pwldog", "dbldog")
	globnames <- c("pure", "geom. line search", "quad. line search", "Powell trust region", "dbl. trust region")
	
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
			cat("Begin", paste(g,l), "\n")
		
		mytime <- system.time(
							  res <- GNE.nseq(xinit, ..., method=method, 
										  control=c(control, list(LM.param=l)), global=g) 
							  , gcFirst = TRUE)
		res <- c(res, method=method, global=g, LM.param=l, time= mytime[3], name=paste(g,l))
		if(echo)
			cat("End\n")
		reslist <- c(reslist, list(res))
	}	
	
	fullnames2 <- c(fullnames2, "LM - adaptive")
	g <- "none"
	l <- "adaptive"
	if(echo)
		cat("Begin", paste(g,l), "\n")
	
	mytime <- system.time(
						  res <- GNE.nseq(xinit, ..., method=method, 
										  control=c(control, list(LM.param=l)), global=g) 
						  , gcFirst = TRUE)
	res <- c(res, method=method, global=g, LM.param=l, time= mytime[3], name=paste(g,l))
	if(echo)
		cat("End\n")
	reslist <- c(reslist, list(res))
	
	
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




bench.GNE.fpeq <- function(xinit, ..., echo=FALSE, control.outer=list(), 
	control.inner=list())
{
	
	problems <- c("NIR", "VIR")
	merits <- c("NI", "VI", "FP")
	methods <- c("pure", "vH", "UR", "RRE", "MPE", "SqRRE", "SqMPE")
#	methods <- "MPE"
	
	fullnames1 <- paste(rep(methods, each=length(problems)*length(merits)),
						rep(problems, each=length(merits)*length(methods)), 
						rep(merits, length(problems)*length(methods)),
						sep="-")
	idname <- NULL	
	reslist <- list()
	
	for(met in methods)
	for(p in problems)
	for(mer in merits)
	{
		if(echo)
			cat("Begin", paste(met, p, mer), "\n")
		incompat <- (p == "NIR" && mer == "VI")||(p == "VIR" && mer == "NI")||(met == "vH" && mer == "FP")
		if(!incompat)
		{
			idname <- c(idname, TRUE)
			if(met != "UR")
			{
				mytime <- system.time(
					res <- GNE.fpeq(xinit, ..., method=met, problem=p, merit=mer, 
									control.outer=control.outer, control.inner=control.inner) 
								  , gcFirst = TRUE)
			}else
			{
				mytime <- system.time(
					res <- GNE.fpeq(xinit, ..., method=met, problem=p, merit=mer, stepfun=decrstep5,
									control.outer=control.outer, control.inner=control.inner) 
									  , gcFirst = TRUE)
			}
			res <- c(res, problem=p, method=met, merit=mer, time= mytime[3], 
					 name=paste(p, mer, met))
			reslist <- c(reslist, list(res))
		}else
			idname <- c(idname, FALSE)		
		if(echo)
			cat("End\n")

	}	
	
	
	
	fpfncall.outer <- sapply(1:length(reslist), function(i) reslist[[i]]$outer.counts[1] )
	merfncall.outer <- sapply(1:length(reslist), function(i) reslist[[i]]$outer.counts[2] )
	gapfncall.inner <- sapply(1:length(reslist), function(i) reslist[[i]]$inner.counts[1] )
	grgapfncall.inner <- sapply(1:length(reslist), function(i) reslist[[i]]$inner.counts[2] )
	x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
	normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
	comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
	checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )	
	codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )
	
	
	data.frame( method= fullnames1[idname], fpfncall.outer, merfncall.outer, 
			   gapfncall.inner, grgapfncall.inner, comptime, normFx, codes, x )
}






bench.GNE.minpb <- function(xinit, ..., echo=FALSE, control.outer=list(), 
	control.inner=list())
{
	
	problems <- c("NIR", "VIR")
	methods <- c("BB", "BFGS", "CG")
	
	fullnames1 <- paste(rep(methods, each=length(problems)),
						rep(problems, length(methods)), 
						sep="-")
	reslist <- list()
	
	for(met in methods)
	for(p in problems)
	{
		if(echo)
		cat("Begin", paste(met, p), "\n")
		mytime <- system.time(
			res <- GNE.minpb(xinit, ..., method=met, problem=p, 
				control.outer=control.outer, control.inner=control.inner) 
							  , gcFirst = TRUE)
		
		res <- c(res, problem=p, method=met, time= mytime[3], name=paste(p, met))
			
		if(met == "BB")
			print(res)
	
		reslist <- c(reslist, list(res))
		if(echo)
			cat("End\n")
		
	}	
	
	
	
	minfncall.outer <- sapply(1:length(reslist), function(i) reslist[[i]]$outer.counts[1] )
	grminfncall.outer <- sapply(1:length(reslist), function(i) reslist[[i]]$outer.counts[2] )
	gapfncall.inner <- sapply(1:length(reslist), function(i) reslist[[i]]$inner.counts[1] )
	grgapfncall.inner <- sapply(1:length(reslist), function(i) reslist[[i]]$inner.counts[2] )
	x <- t(sapply(1:length(reslist), function(i) reslist[[i]]$par ))
	normFx <- sapply(1:length(reslist), function(i) reslist[[i]]$value )
	comptime <- sapply(1:length(reslist), function(i) reslist[[i]]$time )
	checknam <- sapply(1:length(reslist), function(i) reslist[[i]]$name )	
	codes <- sapply(1:length(reslist), function(i) reslist[[i]]$code )
	
	
	print(sapply(1:length(reslist), function(i) reslist[[i]]$outer.counts))
	print(length(fullnames1))
	print(minfncall.outer)
	print(grminfncall.outer)
	print(unlist(minfncall.outer))
	print(unlist(grminfncall.outer))
	print(length(gapfncall.inner))
	print(length(grgapfncall.inner))
	print(length(comptime))
	print(length(normFx))
	print(length(codes))
	print(length(x))
	
	data.frame( method= fullnames1, minfncall.outer, grminfncall.outer, 
			   gapfncall.inner, grgapfncall.inner, comptime, 
			   normFx, codes, x )
}




