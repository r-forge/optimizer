#functions of the Variational Inequality Reformulation of the GNEP


gapVIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	#sanity check		
	arg <- testarggapVIR(x, dimx, grobj, arggrobj)
	
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
	grobjregi <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		xy <- x
		xy[idx_i] <- y[idx_i]
		norm_xy <- sum( (x[idx_i] - y[idx_i])^2 )
		arg$obj(x, i, i, arg$arggrobj) - arg$grobj(xy, i, i, arg$arggrobj) - par$alpha/2*norm_xy
	}
	sum( sapply(1:nplayer, grobjregi) )
}


gradxgapVIR <- function(x, y, dimx, grobj, arggrobj, 
	heobj, argheobj, param=list(), echo=FALSE)
{
	arg <- testarggradgapVIR(x, dimx, grobj, arggrobj, 
		heobj, argheobj, echo)

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
	

	stop("not yet implemented")
	
}

gradygapVIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	arg <- testarggradgapVIR(x, dimx, grobj, arggrobj, echo)
	
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

	stop("not yet implemented")
}

fpVIR <- function(x, dimx, obj, argobj, joint, argjoint,  
	grobj, arggrobj, jacjoint, argjacjoint, param=list(), 
	echo=FALSE, control=list(), yinit=NULL, optim.method="default")
{
	arg <- testargfpVIR(x, dimx, obj, argobj, joint, argjoint,  
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
	
	
		stop("not yet implemented")
		
	
}

