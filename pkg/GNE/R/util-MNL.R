#####################
#price ratio function


PROBj2kRatio <- function(x, j, k, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	
#	if(any(x == 0) || any(x < 0))
#		res <- ifelse(j == k, 0, 1)
#	else if(is.infinite(x[k]) || is.infinite(x[j]))
#		res <- 0
#	else   
	res <- .Call("doPROBj2kRatio", as.double(x), as.integer(j), 
				 as.integer(k), as.double(as.numeric(param_j)))
	res
}



GrPROBj2kRatio <- function(x, j, k, ideriv, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	if(any(is.infinite(x)))
		stop("x has infinite values.")
#	if(any(x == 0) || any(x < 0))
#		res <- ifelse(j == k, 0, 1)
#	else   
	res <- .Call("doGradPROBj2kRatio", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.double(param_j))

	res
}

GrGrPROBj2kRatio <- function(x, j, k, ideriv, mderiv, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	if(any(is.infinite(x)))
		stop("x has infinite values.")
#	if(any(x == 0) || any(x < 0))
#		res <- ifelse(j == k, 0, 1)
#	else   
	res <- .Call("doGradGradPROBj2kRatio", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.integer(mderiv), as.double(param_j))
	
	res
}


#####################
#price diff function


PROBj2kDiff <- function(x, j, k, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	
#	if(any(x == 0) || any(x < 0))
#		res <- ifelse(j == k, 0, 1)
#	else if(is.infinite(x[k]) || is.infinite(x[j]))
#		res <- 0
#	else   
	res <- .Call("doPROBj2kDiff", as.double(x), as.integer(j), 
				 as.integer(k), as.double(as.numeric(param_j)))
	
	res
}



GrPROBj2kDiff <- function(x, j, k, ideriv, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	if(any(is.infinite(x)))
		stop("x has infinite values.")
#	if(any(x == 0) || any(x < 0))
#	res <- ifelse(j == k, 0, 1)
#	else   
	res <- .Call("doGradPROBj2kDiff", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.double(param_j))
	
	res
}

GrGrPROBj2kDiff <- function(x, j, k, ideriv, mderiv, param_j)
{
	if(!is.numeric(x))
		stop("wrong argument x.")
	if(any(is.na(x)))
		stop("x has missing values.")
	if(any(is.infinite(x)))
		stop("x has infinite values.")
#	if(any(x == 0) || any(x < 0))
#	res <- ifelse(j == k, 0, 1)
#	else   
	res <- .Call("doGradGradPROBj2kDiff", as.double(x), as.integer(j), 
				 as.integer(k), as.integer(ideriv), as.integer(mderiv), as.double(param_j))
	
	res
}

