# geometric backtracking line search - minstep based on Dennis and Schnabel
linesearch.geom <- function(xk, dk, slopek, con, merit, ...)
{
	normfk <- merit(xk, ...)
	stepk <- 1
			
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
			
	while(stepk > minstep)
	{
		normfp <- merit(xk + stepk*dk, ...)
		
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
		
		#check Armijo condition
		if(normfp <= normfk + con$btol * stepk * slopek)
		{
			break
		}
		
		inner.iter <- inner.iter + 1	
		stepk <- con$sigma * stepk
	}
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)

}			

linesearch.geom.cond <- function(xk, dk, slopek, con, merit, checkint, ...)
{
	normfk <- merit(xk, ...)
	stepk <- 1
			
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
			
	while(stepk > minstep)
	{
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\n")			
		
		if(checkint(xk + stepk*dk, ...))
		{
			normfp <- merit(xk + stepk*dk, ...)
	
			#traces in R console	
			if(con$trace >= 3)
				cat("left term", normfp, "\tright term\t", 
					normfk + con$btol * stepk * slopek, "descent", stepk * slopek, "\n")			

			#check Armijo condition
			if(normfp <= normfk + con$btol * stepk * slopek)
			{
				break
			}
		}
		
		inner.iter <- inner.iter + 1	
		stepk <- con$sigma * stepk
	}
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)

}			



# quadratic backtracking line search - minstep based on Dennis and Schnabel
linesearch.quad <- function(xk, dk, slopek, con, merit, ...)
{
	normfk <- merit(xk, ...)
	stepk <- 1
			
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
			
	while(stepk > minstep)
	{
		normfp <- merit(xk + stepk*dk, ...)
		
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
		
		
		#check Armijo condition
		if(normfp <= normfk + con$btol * stepk * slopek)
			break
		
		inner.iter <- inner.iter + 1	
		stepk <- - as.numeric( (stepk)^2 * slopek / 2 / (normfp - normfk - slopek)	)	

	}	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)
}			
	
	