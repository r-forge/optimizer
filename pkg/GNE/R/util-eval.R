evalwitharglist <- function(f, x, addarg)
{
	if(length(addarg) == 0)
	{	
		f(x)
	}else
	{	
		do.call(f, c(list(x), addarg) )
	}	
}


catlist <- function(x, ..., issummary=TRUE)
{
	if(issummary)
		cat(capture.output(summary(x)), sep="\n", ...)
	else
		cat(capture.output(print(x)), sep="\n", ...)
}



list2array <- function(x, listname=1:length(x))
array(unlist(x), 
	dim = c(dim(x[[1]]), length(x)),
	dimnames=c(dimnames(x[[1]]), list(listname))
) 

testfunc <- function(f, x=NULL, arg=NULL, echo=TRUE, errmess="")
{
	nbtotarg <- length(formals(f))
	if(is.logical(echo))
		echo <- 1
	if(echo >= 2)
		print(formals(f))
	
	if(!is.null(x))
	{
		if(!is.null(arg))
		{
			if(nbtotarg < 2)
				stop(errmess)
			else if(nbtotarg == 2)
				arglist <- list(x, arg)
			else if(nbtotarg >= 3)
				arglist <- c(list(x), as.list(rep(1, nbtotarg-2)), list(arg))
		}else
		{
			if(nbtotarg < 1)
				stop(errmess)
			else if(nbtotarg == 1)
				arglist <- list(x)
			else if(nbtotarg >= 2)
				arglist <- c(list(x), as.list(rep(1, nbtotarg-1)))
		}
	}else
	{
		if(!is.null(arg))
		{
			if(nbtotarg < 1)
				stop(errmess)
			else if(nbtotarg == 1)
				arglist <- list(arg)
			else if(nbtotarg >= 2)
				arglist <- c(as.list(rep(1, nbtotarg-1)), list(arg))
		}else
		{
			if(nbtotarg == 0)
				arglist <- list()
			else if(nbtotarg >= 1)
				arglist <- as.list(rep(1, nbtotarg))
		}		
	}
	test.try <- try(do.call(f, arglist), silent=echo >= 2)
					
	if(class(test.try) == "try-error")
	{
		if(echo >= 1)
		{
			cat("Error when calling function, below the try output.\n")
			print(test.try)
			cat("Arguments are:\n")
			print(arglist)
			stop(errmess)
		}
	}
	invisible()
}


