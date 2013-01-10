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
		res <- GNE.fpeq(xinit, method=method, control.outer=control, ...)
	
	if(approach == "minimization")
		res <- GNE.min(xinit, method=method, control=control, ...)
	
	c(res, list(approach=approach))
}

GNE.nseq <- function(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	compl, gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="default", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "Newton"
	
	argtest1 <- testargfunSSR(init, dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
						 compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint)
	
	#basic tests for funSSR
	test.try <- try( funSSR(init, dimx, dimlam, grobj, arggrobj, constr, argconstr,  
							grconstr, arggrconstr, compl, argcompl, dimmu, joint, argjoint,
							grjoint, arggrjoint), silent=silent )
	
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(init) has infinite or NaN values.", fvec=NA) )

	argtest2 <- testargjacSSR(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
						 heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint,
						 hejoint, arghejoint)	
	
	#basic tests for jacSSR
	test.try <- try( jacSSR(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
					  heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint,
					  grjoint, arggrjoint, hejoint, arghejoint), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac Phi(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac Phi(init) has infinite or NaN values.", fvec=NA) )
	
	#wrapped functions
	myfunSSR <- function(x, argfun, argjac)
		evalwitharglist(funSSR, x, argfun) 
	
	myjacSSR <- function(x, argfun, argjac)
		evalwitharglist(jacSSR, x, argjac)
	
	arg1 <- list(argtest1$dimx, argtest1$dimlam, argtest1$grobj, argtest1$arggrobj, 
				 argtest1$constr, argtest1$argconstr, argtest1$grconstr, argtest1$arggrconstr, 
				 argtest1$compl, argtest1$argcompl, argtest1$dimmu, 
				 argtest1$joint, argtest1$argjoint, argtest1$grjoint, argtest1$arggrjoint)
	arg2 <- list(argtest2$dimx, argtest2$dimlam, argtest2$heobj, argtest2$argheobj, 
				 argtest2$constr, argtest2$argconstr, argtest2$grconstr, argtest2$arggrconstr, 
				 argtest2$heconstr, argtest2$argheconstr, argtest2$gcompla, argtest2$gcomplb, argtest2$argcompl, 
				 argtest2$dimmu, argtest2$joint, argtest2$argjoint, argtest2$grjoint, 
				 argtest2$arggrjoint, argtest2$hejoint, argtest2$arghejoint)	
	if(!silent)
		print("init completed.")
	
	res <- nseq(init, myfunSSR, myjacSSR, argfun=arg1, argjac=arg2, method=method, 
				control=control, silent=silent, ...)	
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	
	res
}

GNE.ceq <- function(init, dimx, dimlam, grobj, arggrobj, heobj, argheobj, 
	constr, argconstr, grconstr, arggrconstr, heconstr, argheconstr,
	dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint, 
	method="PR", control=list(), silent=TRUE, ...)
{
	if(method == "default") method <- "PR"
	
	argtest1 <- testargfunCER(init, dimx, dimlam, grobj, arggrobj, constr, argconstr, 
							  grconstr, arggrconstr, dimmu, joint, argjoint, grjoint, arggrjoint)
	
	#basic tests for funCER
	test.try <- try( funCER(init, dimx, dimlam, grobj, arggrobj, constr, argconstr, 
				grconstr, arggrconstr, dimmu, joint, argjoint, grjoint, arggrjoint), silent=silent )

	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate H(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="H(init) has infinite or NaN values.", fvec=NA) )
	
	argtest2 <- testargjacCER(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, 
							  arggrconstr, heconstr, argheconstr, dimmu, joint, argjoint, grjoint, 
							  arggrjoint, hejoint, arghejoint)
	
	#basic tests for jacCER
	test.try <- try( jacCER(init, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, 
							arggrconstr, heconstr, argheconstr, dimmu, joint, argjoint,
							grjoint, arggrjoint, hejoint, arghejoint), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac H(init).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac H(init) has infinite or NaN values.", fvec=NA) )
	
	#wrapped functions
	myfunCER <- function(x, argfun, argjac)
		evalwitharglist(funCER, x, argfun) 
	
	myjacCER <- function(x, argfun, argjac)
		evalwitharglist(jacCER, x, argjac)
	
	arg1 <- list(dimx = argtest1$dimx, dimlam = argtest1$dimlam, grobj = argtest1$grobj, 
		arggrobj = argtest1$arggrobj, constr = argtest1$constr, argconstr = argtest1$argconstr, 
		grconstr = argtest1$grconstr, arggrconstr = argtest1$arggrconstr, dimmu = argtest1$dimmu, 
		joint = argtest1$joint, argjoint = argtest1$argjoint, grjoint = argtest1$grjoint, 
		arggrjoint = argtest1$arggrjoint)
	arg2 <- list(argtest2$dimx, argtest2$dimlam, argtest2$heobj, argtest2$argheobj, 
				 argtest2$constr, argtest2$argconstr, argtest2$grconstr, argtest2$arggrconstr, 
				 argtest2$heconstr, argtest2$argheconstr, argtest2$dimmu, argtest2$joint, 
				 argtest2$argjoint, argtest2$grjoint, argtest2$arggrjoint,
				 argtest2$hejoint, argtest2$arghejoint)	
	if(!silent)
		print("init completed.")
		
	res <- ceq(init, dimx, dimlam, Hfinal=myfunCER, jacHfinal=myjacCER, 
			   argfun=arg1, argjac=arg2, method=method, control=control, 
			   silent=silent, ...)	
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	
	res
}


GNE.fpeq <- function(init, dimx, obj, argobj, grobj, arggrobj, 
	heobj, argheobj, joint, argjoint, jacjoint, argjacjoint, 
	method = "default", problem = c("NIR", "VIR"), 
	merit = c("NI", "VI", "FP"), order=1, control.outer=list(), 
	control.inner=list(), silent=TRUE, param=list(), stepfunc, argstep, ...)
{

	if(method == "default") method <- "MPE"
	method <- match.arg(method, c("pure", "vH", "UR", "RRE", "MPE", "SqRRE", "SqMPE"))
	problem <- match.arg(problem, c("NIR", "VIR"))
	merit <- match.arg(merit, c("NI", "VI", "FP"))
	
	if(method %in% c("vH") && merit == "FP")
		stop("incompatible method and merit arguments.")
	
	if(order > 3 || order < 1)
		stop("wrong order argument.")
	if(problem == "NIR")
	{
		arggap <- testarggapNIR(init, dimx, obj, argobj)
		argfp <- testargfpNIR(init, dimx, obj, argobj, joint, argjoint,  
					 grobj, arggrobj, jacjoint, argjacjoint)
		
		yfun <- function(x)
			fpNIR(x, argfp$dimx, argfp$obj, argfp$argobj, argfp$joint, 
				  argfp$argjoint, argfp$grobj, argfp$arggrobj, argfp$jacjoint,
				  argfp$argjacjoint, echo=!silent, control=control.inner, 
				  param=param, ...)

		if(merit == "FP")
			merit <- NULL
		else if(merit == "NI")
		{
			merit <- function(x, y=NULL)
			{
				if(is.null(y))
				{
					y <- yfun(x)
					res <- gapNIR(x, y$par, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent)
					y$iter <- ifelse(!is.na(y$counts[2]), y$counts[2], y$counts[1])
					list(value=res, counts=y$counts+c(1,0), iter=y$iter)
				}else
					list(value=gapNIR(x, y, arggap$dimx, arggap$obj, arggap$argobj, 
								  param=param, echo=!silent),
						 counts=c(1,0), iter=0)
			}
		
		}else
			stop("wrong merit function")
	}
		
	if(!silent)
		print("init completed.")
	res <- fpeq(init, yfun, merit, method, control=control.outer, stepfunc=stepfunc,
				argstep=argstep, silent=silent, order=order, ...)
	class(res) <- "GNE"
	if(!silent)
		print("computation completed.")
	res
}


GNE.min <- function(init, gap, gradgap, arggap=list(), arggrad=list(), method, 
	control=list(), ...)
{
	if(method == "default")
		method <- "BB"
	
	stop("not yet implemented")
	
	res <- minpb(init, gap, gradgap, arggap, arggrad, method, control, ...)
	class(res) <- "GNE"
	res
}





#print function
print.GNE <- function(x, ...)
{
	if (!inherits(x, "GNE"))
		stop("Use only with 'GNE' objects")
	cat("GNE:", x$par, "\nwith optimal norm", x$value, "\n")
	cat("after ", x$iter, "iterations with exit code", x$code, ".\n")
	cat("Output message:", x$message, "\n")
	if(!is.null(x$counts))	
		cat("Function/grad/hessian calls:", x$counts, "\n")
	if(!is.null(x$outer.counts))	
		cat("Function/grad/hessian calls:", x$outer.counts, x$inner.counts, "\n")	
	if(!is.null(x$fvec))
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

