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
	


	test.try <- try( grobj(z, 1, 1, arggrobj) , silent=FALSE)
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to grobj(z, 1, 1, arggrobj) does not work.", "arguments are", 
					 paste(names(formals(grobj)), collapse=","), ".")
		stop(str)
	}
	test.try <- try( constr(z, 1, argconstr) , silent=FALSE)
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to constr(z, 1, argconstr) does not work.", "arguments are", 
					 paste(names(formals(constr)), collapse=","), ".")
		stop(str)
	}
	test.try <- try( grconstr(z, 1, 1, arggrconstr) , silent=FALSE)
	if(class(test.try) == "try-error")
	{
		str <- paste("the call to grconstr(z, 1, 1, arggrconstr) does not work.", "arguments are", 
					 paste(names(formals(grconstr)), collapse=","), ".")
		stop(str)
	}
	
	
	res <- .Call("doPhi", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 grobj, arggrobj, 
				 constr, argconstr, 
				 grconstr, arggrconstr, 
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
		stop("wrong dimension 1")
	if(length(z) != sum(dimx) + sum(dimlam))
		stop("wrong dimension 3")
	

	
	res <- .Call("doJacPhi", nplayer, z, as.integer(dimx), as.integer(dimlam), 
				 heobj, argheobj, 
				 constr, argconstr, 
				 grconstr, arggrconstr, 
				 heconstr, argheconstr,
				 gcompla, gcomplb, new.env())

	
	res
}


