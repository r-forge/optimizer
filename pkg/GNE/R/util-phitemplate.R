Phi <- function(z, dimx, dimlam,
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	compl)
{
	#sanity check	
	nplayer <- length(dimx)
	if(nplayer != length(dimlam) || length(z) != sum(dimx) + sum(dimlam))
		stop("incompatible dimension for dimlam and dimx")
	if(!is.function(grobj) || !is.function(constr) || !is.function(grconstr) || !is.function(compl))
		stop("one (at least) of these arguments is not a function grobj, constr, grconstr, compl")
	
	#objective gradient
	if(!missing(arggrobj))
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
	if(!missing(argconstr))
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
	if(!missing(arggrconstr))
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
	
	
	res <- .Call("doPhi", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 grobjfinal, arggrobj, 
				 constrfinal, argconstr, 
				 grconstrfinal, arggrconstr, 
				 compl, new.env())

	res
}


JacPhi <- function(z, dimx, dimlam,
	heobj, argheobj, 
	constr, argconstr, 
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	gcompla, gcomplb)
{
	
	nplayer <- length(dimx)
	if(nplayer != length(dimlam))
		stop("wrong dimension.")
	if(length(z) != sum(dimx) + sum(dimlam))
		stop("wrong dimension")
	
	
	#objective hessian
	if(!missing(argheobj))
	{
		test.try <- try( heobj(z, 1, 1, 1, argheobj) , silent=FALSE)
		strcall <- "heobj(z, 1, 1, 1, arggrobj)"
		heobjfinal <- heobj
	}else
	{
		test.try <- try( heobj(z, 1, 1, 1) , silent=FALSE)
		heobjfinal <- function(z, i, j, k, ...)
			heobj(z, i, j, k)
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
	if(!missing(argconstr))
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
	if(!missing(arggrconstr))
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
	
	#constraint hessian
	if(!missing(argheconstr))
	{
		test.try <- try( heconstr(z, 1, 1, 1, argheconstr) , silent=FALSE)
		strcall <- "heconstr(z, 1, 1, 1, argheconstr)"
		heconstrfinal <- heconstr
	}else
	{
		test.try <- try( heconstr(z, 1, 1, 1) , silent=FALSE)
		heconstrfinal <- function(z, i, j, k, ...)
			heconstr(z, i, j, k)
		argheconstr <- list()				   
		strcall <- "heconstr(z, 1, 1, 1)"	
	}	
	if(class(test.try) == "try-error")
	{
		print(argheconstr)
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(argheconstr)), collapse=","), ".")
		stop(str)
	}
	

	
	res <- .Call("doJacPhi", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 heobjfinal, argheobj, 
				 constrfinal, argconstr, 
				 grconstrfinal, arggrconstr, 
				 heconstrfinal, argheconstr,
				 gcompla, gcomplb, new.env())

	
	res
}


