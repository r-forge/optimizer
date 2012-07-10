#functions of the SemiSmooth Reformulation of the GNEP
#z = (x , lambda)



#function phi of the SSR
funSSR <- function(z, dimx, dimlam,
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	compl, argcompl)
{
	#sanity check	
	nplayer <- length(dimx)
	if(nplayer != length(dimlam) || length(z) != sum(dimx) + sum(dimlam))
		stop("incompatible dimension for dimlam and dimx")
	if(!is.function(grobj) || !is.function(constr) || !is.function(grconstr) || !is.function(compl))
		stop("one (at least) of these arguments is not a function grobj, constr, grconstr, compl")
	
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
	
	#complementarity function
	if(!missing(argcompl) && !is.null(argcompl))
	{
		complfinal <- function(a, b, argcompl)
			evalwitharglist(compl, a, c(list(b), argcompl))
		
		test.try <- try( complfinal(1, 1, argcompl) , silent=FALSE)
		strcall <- "complfinal(1, 1, argcompl)"
	}else
	{
		test.try <- try( compl(1, 1) , silent=FALSE)
		complfinal <- function(a, b, argcompl)
			compl(a, b)
		argcompl <- list()
		strcall <- "compl(1, 1)"
	}
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall, "does not work.", "arguments are", 
					 paste(names(formals(compl)), collapse=","), ".")
		stop(str)
	}
	
	res <- .Call("dofunSSR", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 grobjfinal, arggrobj, constrfinal, argconstr, 
				 grconstrfinal, arggrconstr, complfinal, argcompl, new.env())

	res
}


#Jacobian of phi of the SSR
jacSSR <- function(z, dimx, dimlam,
	heobj, argheobj, 
	constr, argconstr, 
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	gcompla, gcomplb, argcompl)
{
	
	nplayer <- length(dimx)
	if(nplayer != length(dimlam))
		stop("wrong dimension.")
	if(length(z) != sum(dimx) + sum(dimlam))
		stop("wrong dimension")
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
	
	#constraint hessian
	if(!missing(argheconstr) && !is.null(argheconstr))
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
		str <- paste("the call to", strcall,"does not work.", "arguments are", 
					 paste(names(formals(argheconstr)), collapse=","), ".")
		stop(str)
	}
	
	#complementarity function
	if(!missing(argcompl) && !is.null(argcompl))
	{
		gcomplafinal <- function(a, b, argcompl)
			evalwitharglist(gcompla, a, c(list(b), argcompl))
		gcomplbfinal <- function(a, b, argcompl)
			evalwitharglist(gcomplb, a, c(list(b), argcompl))
		
		test.try <- try( gcomplafinal(1, 1, argcompl) , silent=FALSE)
		strcall <- "gcomplafinal(1, 1, argcompl)"
		
	}else
	{
		test.try <- try( gcompla(1, 1) , silent=FALSE)
		gcomplafinal <- function(a, b, argcompl)
			gcompla(a, b)	
		gcomplbfinal <- function(a, b, argcompl)
			gcomplb(a, b)	
		argcompl <- list()
		strcall <- "gcompla(1, 1)"
	}
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to", strcall, "does not work.", "arguments are", 
					 paste(names(formals(gcompla)), collapse=","), ".")
		stop(str)
	}
	

	
	res <- .Call("dojacSSR", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 heobjfinal, argheobj, constrfinal, argconstr, 
				 grconstrfinal, arggrconstr, heconstrfinal, argheconstr,
				 gcomplafinal, gcomplbfinal, argcompl, new.env())

	
	res
}


