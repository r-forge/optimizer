#functions of the Nikaido Isoda Reformulation of the GNEP


gapNIR <- function(x, y, dimx, obj, argobj, param=list(), echo=FALSE)
{
	#sanity check		
	arg <- testarggapNIR(x, dimx, obj, argobj)
	
	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")

	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	
	#ith increment
	objregi <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		xy <- x
		xy[idx_i] <- y[idx_i]
		norm_xy <- sum( (x[idx_i] - y[idx_i])^2 )
		arg$obj(x, i, arg$argobj) - arg$obj(xy, i, arg$argobj) - par$alpha/2*norm_xy
	}
	sum( sapply(1:nplayer, objregi) )
}


gradxgapNIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	arg <- testarggradgapNIR(x, dimx, grobj, arggrobj, echo)

	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")
	
	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	
	grad_obji <- function(x, i)
		sapply(1:sum(dimx), function(j) arg$grobj(x, i, j, arg$arggrobj))
	deltagrad <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		xy <- x
		xy[idx_i] <- y[idx_i]
		grad_obji(x, i) - grad_obji(xy, i)
	}
	gradi_obji <- function(i)	
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		xy <- x
		xy[idx_i] <- y[idx_i]
		sapply(idx_i, function(j) arg$grobj(xy, i, j, arg$arggrobj))
	}
	
	res <- rowSums(sapply(1:nplayer, deltagrad)) 
	res <- res + unlist(sapply(1:nplayer, gradi_obji))
	res - par$alpha*(x-y)
}

gradygapNIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	arg <- testarggradgapNIR(x, dimx, grobj, arggrobj, echo)
	
	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")
	
	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	
	gradi_obji <- function(i)	
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		xy <- x
		xy[idx_i] <- y[idx_i]
		sapply(idx_i, function(j) arg$grobj(xy, i, j, arg$arggrobj))
	}
	- unlist(sapply(1:nplayer, gradi_obji)) - par$alpha*(x-y)
}

fpNIR <- function(x, dimx, obj, argobj, joint, argjoint,  
	grobj, arggrobj, jacjoint, argjacjoint, param=list(), 
	echo=FALSE, control=list(), yinit=NULL, optim.method="default")
{
	arg <- testargfpNIR(x, dimx, obj, argobj, joint, argjoint,  
			grobj, arggrobj, jacjoint, argjacjoint, echo)

	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	if(optim.method == "default")
		optim.method <- "BFGS"
	#default control parameters
	con1 <- list(trace = echo, eps=1e-6, method=optim.method)
	namc1 <- names(con1)
	con1[namc1 <- names(control)] <- control
	con2 <- list(trace = echo, reltol=1e-6)
	namc2 <- names(con2)
	con2[namc2 <- names(control)] <- control
	
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(is.null(yinit) && !is.null(arg$joint))
	{
		yinit <- rejection(function(x) arg$joint(x, arg$argjoint), n, method="norm")
		if(echo)
			cat("initial point in fpNIR", yinit, "-", arg$joint(yinit, arg$argjoint), "\n")
	}else if(is.null(yinit) && is.null(arg$joint))
	{
		yinit <- rnorm(n)
	}else
	{
		yinit <- yinit[1:n]
	}
	
	
	if(is.null(arg$joint) && !is.null(arg$grobj))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapNIR(z, y, dimx, arg$obj, arg$argobj, param)
			gr <- function(y, z, param) -gradygapNIR(z, y, dimx, arg$grobj, arg$arggrobj, param) 
			
			res <- optim(yinit, fn=fn, gr=gr, method=optim.method, 
								  z=x, param=par, control=con2)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.")
			c(res, optim.function="optim", optim.method=optim.method)
		}
	}else if(is.null(arg$joint) && is.null(arg$grobj))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapNIR(z, y, dimx, arg$obj, arg$argobj, param)
			
			res <- optim(yinit, fn=fn, method=optim.method, 
						 z=x, param=par, control=con2)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.")
			c(res, optim.function="optim", optim.method=optim.method)
		}
	}else if(!is.null(arg$joint) && !is.null(arg$jacjoint))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapNIR(z, y, dimx, arg$obj, arg$argobj, param)
			gr <- function(y, z, param) -gradygapNIR(z, y, dimx, arg$grobj, arg$arggrobj, param) 
			hin <- function(y, z, param) -arg$joint(y, arg$argjoint)
			hin.jac <- function(y, z, param) -arg$jacjoint(y, arg$argjacjoint)
			
			res <- constrOptim.nl(yinit, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac, 
						   z=x, param=par, control.outer=con1)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.")
			c(res, optim.function="constrOptim.nl", optim.method=optim.method)
		}
	}else if(!is.null(arg$joint) && is.null(arg$jacjoint))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapNIR(z, y, dimx, arg$obj, arg$argobj, param)
			hin <- function(y, z, param) -arg$joint(y, arg$argjoint)
			
			res <- constrOptim.nl(yinit, fn=fn, hin=hin, 
								  z=x, param=par, control.outer=con1)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.")
			c(res, optim.function="constrOptim.nl", optim.method=optim.method)
		}
	}else
		stop("missing argument in obj, grobj, joint, jacjoint.")
	yofx(x)
}

