#rejection algorithm with different distributions

rejection <- function(constr, nvars, LB=0, UB=1, ..., echo=FALSE, 
	method=c("unif","norm", "normcap"), control=list())
{
	#default control parameters
	con <- list(mean=(LB + UB)/2, sd=-(UB-LB) /4 /qnorm(0.025))
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	
	#adjust the lentgh of lower/upper bounds
	LB <- rep(LB, length=nvars)	
	method <- match.arg(method, c("unif","norm","normcap"))
	
	if(method == "unif")
	{
		#draw
		x <- runif(nvars, LB, UB)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) )
		{
			if(echo) print(x)
			x <- runif(nvars, LB, UB)
		}
		return(x)
	}
	
	if(method == "norm")
	{
		#draw
		x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) )
		{
			if(echo) print(x)
			x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		}
		return(x)
	}
	
	if(method == "normcap")
	{
		#draw
		x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) || sum( (UB-x)^2 ) <= 1)
		{
			if(echo) print(x)
			x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		}
		return(x)
	}
	
	stop("wrong method for rejection algo.")
	
}

	