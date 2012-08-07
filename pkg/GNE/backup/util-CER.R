#functions of the Constrained Equation Reformulation of the GNEP
#z = (x, lam, w)
#with size (n, m, m)


funCER <- function(z, dimx, dimlam, dimw,
grobj, arggrobj, 
constr, argconstr,  
grconstr, arggrconstr)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	nplayer <- length(dimx)
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]

	#sanity check	
	if(nplayer != length(dimlam) || length(z) != sum(dimx) + sum(dimlam) + sum(dimw) || sum(dimlam) != sum(dimw))
		stop("incompatible dimension for dimlam, dimx, dimw")
	if(!is.function(grobj) || !is.function(constr) || !is.function(grconstr))
		stop("one (at least) of these arguments is not a function grobj, constr, grconstr")
	
	#objective gradient
	if(!missing(arggrobj) && !is.null(arggrobj))
	{
		test.try <- try( grobj(z, 1, 1, arggrobj) , silent=FALSE)
		strcall <- "grobj(z, 1, 1, arggrobj)"
		grobjfinal <- grobj
	}else
	{
		test.try <- try( grobj(z, 1, 1) , silent=FALSE)
		grobjfinal <- function(z, i, j, ...)
			grobj(z, i, j)
		arggrobj <- list()				   
		strcall <- "grobj(z, 1, 1)"	
	}	
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(grobj)), collapse=","), ".")
		stop(str)
	}

#constraint	
	if(!missing(argconstr) && !is.null(argconstr))
	{
		test.try <- try( constr(z, 1, argconstr) , silent=FALSE)
		strcall <- "constr(z, 1, argconstr)"
		constrfinal <- constr
	}else
	{
		test.try <- try( constr(z, 1) , silent=FALSE)
		constrfinal <- function(z, i, ...)
			constr(z, i)
		argconstr <- list()	
		strcall <- "constr(z, 1)"
	}
	
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(constr)), collapse=","), ".")
		stop(str)
	}
	
#constraint gradient
	if(!missing(arggrconstr) && !is.null(arggrconstr))
	{
		test.try <- try( grconstr(z, 1, 1, arggrconstr) , silent=FALSE)
		strcall <- "grconstr(z, 1, 1, arggrconstr)"
		grconstrfinal <- grconstr
	}else
	{
		test.try <- try( grconstr(z, 1, 1) , silent=FALSE)
		grconstrfinal <- function(z, i, j, ...)
			grconstr(z, i, j)
		arggrconstr <- list()
		strcall <- "grconstr(z, 1, 1)"
	}
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall, "does not work.", "arguments are", 
					 paste(names(formals(grconstr)), collapse=","), ".")
		stop(str)
	}
	
	
#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) )
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )

	GrLagri <- function(i) 
	{
#i index for player, j index for x_ij
		sapply(index4x[1,i]:index4x[2,i], function(j)			   
			   grobjfinal(x, i, j, arggrobj) + lam[index4lam[1,i]:index4lam[2,i]] %*% grconstrfinal(x, i, j, arggrobj)
			   ) 
	}
	Constri <- function(i)
		constrfinal(z, i, argconstr)
	
	c(	unlist(sapply(1:nplayer, GrLagri)),
		unlist(sapply(1:nplayer, Constri)) + w,
		lam * w )
}



#z = (x, lam, w)
#with size (n, m, m)
jacCER <- function(z, dimx, dimlam, dimw,
heobj, argheobj, 
constr, argconstr, 
grconstr, arggrconstr, 
heconstr, argheconstr
)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	nplayer <- length(dimx)
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[-(1:(n+m))]
	
	#sanity check	
	nplayer <- length(dimx)
	if(nplayer != length(dimlam) || length(z) != sum(dimx) + sum(dimlam) + sum(dimw))
		stop("incompatible dimension for dimlam, dimx, dimw")
	if(!is.function(heobj) || !is.function(constr) || !is.function(grconstr) || !is.function(heconstr))
		stop("one (at least) of these arguments is not a function heobj, constr, grconstr, heconstr")
	
	
	
#objective hessian
	if(!missing(argheobj) && !is.null(argheobj))
	{
		test.try <- try( heobj(z, 1, 1, 1, argheobj) , silent=FALSE)
		strcall <- "heobj(z, 1, 1, 1, arggrobj)"
		heobjfinal <- heobj
	}else
	{
		test.try <- try( heobj(z, 1, 1, 1) , silent=FALSE)
		heobjfinal <- function(z, i, j, k, ...) heobj(z, i, j, k)
		argheobj <- list()				   
		strcall <- "heobj(z, 1, 1, 1)"	
	}	
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(argheobj)), collapse=","), ".")
		stop(str)
	}
	
#constraint	
	if(!missing(argconstr) && !is.null(argconstr))
	{
		test.try <- try( constr(z, 1, argconstr) , silent=FALSE)
		strcall <- "constr(z, 1, argconstr)"
		constrfinal <- constr
	}else
	{
		test.try <- try( constr(z, 1) , silent=FALSE)
		constrfinal <- function(z, i, ...) constr(z, i)
		argconstr <- list()	
		strcall <- "constr(z, 1)"
	}
	
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(constr)), collapse=","), ".")
		stop(str)
	}
	
#constraint gradient
	if(!missing(arggrconstr) && !is.null(arggrconstr))
	{
		test.try <- try( grconstr(z, 1, 1, arggrconstr) , silent=FALSE)
		strcall <- "grconstr(z, 1, 1, arggrconstr)"
		grconstrfinal <- grconstr
	}else
	{
		test.try <- try( grconstr(z, 1, 1) , silent=FALSE)
		grconstrfinal <- function(z, i, j, ...) grconstr(z, i, j)
		arggrconstr <- list()
		strcall <- "grconstr(z, 1, 1)"
	}
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall, "does not work.", "arguments are", 
					 paste(names(formals(grconstr)), collapse=","), ".")
		stop(str)
	}	
	
#constraint hessian
	if(!missing(argheconstr) && !is.null(argheconstr))
	{
		test.try <- try( heconstr(z, 1, 1, 1, argheconstr) , silent=FALSE)
		strcall <- "heconstr(z, 1, 1, 1, argheconstr)"
		heconstrfinal <- heconstr
	}else
	{
		test.try <- try( heconstr(z, 1, 1, 1) , silent=FALSE)
		heconstrfinal <- function(z, i, j, k, ...) heconstr(z, i, j, k)
		argheconstr <- list()				   
		strcall <- "heconstr(z, 1, 1, 1)"	
	}	
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(argheconstr)), collapse=","), ".")
		stop(str)
	}

#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) )
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )	
	
	
	
	GrjGriLagri <- function(i, j)
	{
#i index for player, j index for x_j, k for x_i_k		
		sapply(index4x[1,i]:index4x[2,i], function(k) 
			   heobjfinal(z, i, j, k, argheobj) + lam[index4lam[1,i]:index4lam[2,i]] %*% heconstrfinal(z, i, j, k, argheconstr)
			   )
	}
	GrjConstri <- function(i) 
	{
#i index for player
		sapply(index4x[1,i]:index4x[2,i], function(j) grconstrfinal(z, i, j, arggrconstr))
	}
	jacconstrij <- function(i, j)
	{
#i index for player, j index for x_j		
		grconstrfinal(z, i, j, arggrconstr)
		
	}
	

	
#Hessian matrix of the Lagrangian
	ggL <- matrix(0, n, n)
	for(i in 1:nplayer)
		ggL[index4x[1,i]:index4x[2,i] , ] <- sapply(1:sum(dimx), function(j) GrjGriLagri(i,j))	
	
#gradient of the constraint function
	gG <- matrix(0, n, m)
	for(i in 1:nplayer)
		gG[index4x[1,i]:index4x[2,i] , index4lam[1,i]:index4lam[2,i]] <- t( GrjConstri(i) )
	
#Jacobian of the constraint functin
	jacG <- matrix(0, m, n)
	for(i in 1:nplayer)
	for(j in 1:sum(dimx))
		jacG[index4lam[1,i]:index4lam[2,i] , j] <- jacconstrij(i,j) 	
	
	
	rbind(cbind(ggL, gG, matrix(0, n, m)),
		  cbind(jacG, diag(m)*0, diag(m)),
		  cbind(matrix(0, m, n), diag(w), diag(lam)) )

	
}


