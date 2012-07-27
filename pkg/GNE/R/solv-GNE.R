GNE <- function(approach = 
	c("non smooth", "fixed point", "minimization", "constrained equation"), 
	method = "default", xinit, control=list(), ...)
{
	approach <- match.arg(approach, c("non smooth", "fixed point", "minimization", "constrained equation"))
	
	if(approach == "non smooth")
		res <- GNE.nseq(xinit, method=method, control=control, ...)

	if(approach == "constr equation")
		res <- GNE.ceq(xinit, method=method, control=control, ...)
	
	if(approach == "fixed point")
		res <- GNE.fpeq(xinit, method=method, control=control, ...)
	
	if(approach == "minimization")
		res <- GNE.min(xinit, method=method, control=control, ...)
	
	c(res, list(approach=approach))
}

GNE.nseq <- function(xinit, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr, 
compl, gcompla, gcomplb, argcompl, method="Newton", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "Newton"
	if(missing(arggrobj)) arggrobj <- NULL
	if(missing(argheobj)) argheobj <- NULL
	if(missing(argconstr)) argconstr <- NULL
	if(missing(arggrconstr)) arggrconstr <- NULL
	if(missing(argheconstr)) argheconstr <- NULL
	if(missing(argcompl)) argcompl <- NULL
	
	#basic tests for funSSR
	test.try <- try( funSSR(xinit, dimx, dimlam, grobj, arggrobj, 
							constr, argconstr, grconstr, arggrconstr, 
							compl, argcompl), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate Phi(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(xinit) has infinite or NaN values.", fvec=NA) )
	
	#basic tests for jacSSR
	test.try <- try( jacSSR(xinit, dimx, dimlam, heobj, argheobj, constr, argconstr, 
							grconstr, arggrconstr, heconstr, argheconstr, gcompla, 
							gcomplb, argcompl), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac Phi(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac Phi(xinit) has infinite or NaN values.", fvec=NA) )
	
	#wrapped functions
	phifinal <- function(x, argPhi, argjac)
		evalwitharglist(funSSR, x, argPhi) 
	
	jacphifinal <- function(x, argPhi, argjac)
		evalwitharglist(jacSSR, x, argjac)
	
	arg1 <- list(dimx, dimlam, grobj, arggrobj, constr, argconstr, grconstr, 
				 arggrconstr, compl, argcompl)
	arg2 <- list(dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, 
				 arggrconstr, heconstr, argheconstr, gcompla, gcomplb, argcompl)	
		
	res <- nseq(xinit, phifinal, jacphifinal, argPhi=arg1, argjac=arg2, method, control, ...)	
	class(res) <- "GNE"
	res
}

GNE.ceq <- function(xinit, dimx, dimlam, dimw, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	method="IP", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "IP"
	if(missing(arggrobj)) arggrobj <- NULL
	if(missing(argheobj)) argheobj <- NULL
	if(missing(argconstr)) argconstr <- NULL
	if(missing(arggrconstr)) arggrconstr <- NULL
	if(missing(argheconstr)) argheconstr <- NULL

#basic tests for funCER
	test.try <- try( funCER(xinit, dimx, dimlam, dimw, grobj, arggrobj, 
							constr, argconstr, grconstr, arggrconstr), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate H(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="H(xinit) has infinite or NaN values.", fvec=NA) )
	
#basic tests for jacCER
	test.try <- try( jacCER(xinit, dimx, dimlam, dimw, heobj, argheobj,  
							constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac H(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac H(xinit) has infinite or NaN values.", fvec=NA) )
	
#wrapped functions
	Hfinal <- function(x, argH, argjac)
		evalwitharglist(funCER, x, argH) 
	
	jacHfinal <- function(x, argH, argjac)
		evalwitharglist(jacCER, x, argjac)
	
	arg1 <- list(dimx, dimlam, dimw, grobj, arggrobj, constr, argconstr, grconstr, arggrconstr)
	arg2 <- list(dimx, dimlam, dimw, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr)	
	
		
	res <- ceq(xinit, dimx, dimlam, dimw, Hfinal, jacHfinal, argH=arg1, argjac=arg2,
		method, control, global="gline", ...)	
	class(res) <- "GNE"
	res
}


GNE.fpeq <- function(xinit, yfunc, func2, argy=list(), argfunc2=list(), method, control=list(), ...)
{
	stop("not yet implemented.")
	if(method == "default")
		method <- "MPE"
		
	res <- fpeq(xinit, yfunc, func2, argy, argfunc2, method, control, ...)
	class(res) <- "GNE"
	res
}


GNE.min <- function(xinit, gap, gradgap, arggap=list(), arggrad=list(), method, control=list(), ...)
{
	if(method == "default")
		method <- "BB"
	
	res <- minpb(xinit, gap, gradgap, arggap, arggrad, method, control, ...)
	class(res) <- "GNE"
	res
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
	
	print(c(control, list(LM.param="adaptive")))
	
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


#print function
print.GNE <- function(x, ...)
{
	if (!inherits(x, "GNE"))
		stop("Use only with 'GNE' objects")
	cat("GNE:", x$par, "\nwith optimal norm", x$value, "\n")
	cat("after ", x$iter, "iterations with exit code", x$code, ".\n")
	cat("Output message:", x$message, "\n")
	cat("Function/grad/hessian calls:", x$counts, "\n")
	cat("Optimal (vector) value:", x$fvec, "\n")
}


#summary function
summary.GNE <- function(object, ...)
{
	structure(object, class = c("summary.GNE", class(object)))	
}


#print function
print.summary.GNE <- function(x, ...)
{	
	if (!inherits(x, "GNE"))
		stop("Use only with 'GNE' objects")
	cat("GNE:", x$par, "\nwith optimal norm", x$value, "\n")
	cat("after ", x$iter, "iterations with exit code", x$code, ".\n")
	if(!is.null(x$playerdead))
		cat("dead players? ", x$playerdead, ".\n")
}

