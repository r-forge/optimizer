PROBj2k <- function(x, j, k, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
	   stop("x has missing values.")
	
	if(any(x == 0) || any(x < 0))
		res <- ifelse(j == k, 0, 1)
	else if(is.infinite(x[k]) || is.infinite(x[j]))
		res <- 0
	else   
		res <- .Call("doPROBj2k", as.double(x), as.integer(j), 
				 as.integer(k), as.double(as.numeric(param_j)))
	   
#	if(any(x < 0))
#		stop("x has negative values.")
	   
	   
	res
}



GrPROBj2k <- function(x, j, k, ideriv, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	if(any(is.infinite(x)))
		stop("x has infinite values.")
	if(any(x == 0) || any(x < 0))
		res <- ifelse(j == k, 0, 1)
	else   
		res <- .Call("doGradPROBj2k", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.double(param_j))

	res
}

GrGrPROBj2k <- function(x, j, k, ideriv, mderiv, param_j)
{
	if(!is.numeric(x))
	stop("wrong argument x.")
	if(any(is.na(x)))
	stop("x has missing values.")
	if(any(is.infinite(x)))
	stop("x has infinite values.")
	if(any(x == 0) || any(x < 0))
	res <- ifelse(j == k, 0, 1)
	else   
	res <- .Call("doGradGradPROBj2k", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.integer(mderiv), as.double(param_j))
	
	res
}

