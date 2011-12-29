require(GNE)

#---------------------------------------------------------------------------------------------------------
# R implementation
#---------------------------------------------------------------------------------------------------------



#______________
# premium ratio based functions, lapse, grlapse, grgrlapse

#examples of parameters c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)
probj2kRatio <- function(x, j, k, paramj)
{
	
	n <- length(x)
	
	denom <- 1 + sum( exp(paramj["mu"] + paramj["alpha"] * x[j] / x[-j]) )
	if(j == k)
	{
		num <- 1
		rtildemin <- 1 - paramj["rmax"] 
	}else
	{
		num <- exp(paramj["mu"] + paramj["alpha"] * x[j] / x[k])
		rtildemin <- paramj["rmin"] / (n-1)
	}
	as.numeric( rtildemin + (paramj["rmax"] - paramj["rmin"]) * num / denom )
}

grprobj2kRatio <- function(x, j, k, ideriv, paramj)
{
	n <- length(x)
	clgjk <- probj2kRatio(x, j, k, paramj)

	if(ideriv == j)
	{
		res <- (-paramj["alpha"]) * sum( sapply( (1:n)[-j], function(l) probj2kRatio(x, j, l, paramj) / x[l] ) ) * clgjk
	}else
	{
		res <- paramj["alpha"] * x[j] / x[ideriv]^2 * probj2kRatio(x, j, ideriv, paramj) * clgjk
	}
	
	if(j != k)
	{
		res <- res + (ideriv == j) * paramj["alpha"] / x[k] * clgjk - (ideriv == k) * paramj["alpha"] * x[j] / x[k]^2 * clgjk
	}
	
	as.numeric((paramj["rmax"] - paramj["rmin"]) * res)
}

grgrprobj2kRatio <- function(x, j, k, ideriv, mderiv, paramj)
{
	n <- length(x)
	clgji <- probj2kRatio(x, j, ideriv, paramj)
	clgjk <- probj2kRatio(x, j, k, paramj)
	grclgji <- grprobj2kRatio(x, j, ideriv, mderiv, paramj)
	grclgjk <- grprobj2kRatio(x, j, k, mderiv, paramj)
	
	if(ideriv == j)
	{
		sum1 <- sum( sapply( (1:n)[-j], function(l) grprobj2kRatio(x, j, l, mderiv, paramj) / x[l] ) )
		sum2 <- sum( sapply( (1:n)[-j], function(l) probj2kRatio(x, j, l, paramj) / x[l] ) )
		
		res <- -paramj["alpha"] * sum1 * clgjk - paramj["alpha"] * sum2 * grclgjk
		res <- res + paramj["alpha"] / x[mderiv]^2 * probj2kRatio(x, j, mderiv, paramj) * clgjk * (j != mderiv)
			
	}else
	{
		res <- paramj["alpha"] * x[j] / x[ideriv]^2 * grclgji * clgjk + paramj["alpha"] * x[j] / x[ideriv]^2 * clgji * grclgjk
		res <- res + paramj["alpha"] / x[ideriv]^2 * clgji * clgjk * (j == mderiv)
		res <- res - 2*paramj["alpha"] * x[j] / x[ideriv]^3 * clgji * clgjk * (ideriv == mderiv)
	}	
	
	if(j != k && ideriv == j)
		res <- res + paramj["alpha"] / x[k] * grclgjk - paramj["alpha"] / x[k]^2 * clgjk * (k == mderiv)
	if(j != k && ideriv == k)
	{	
		res <- res - paramj["alpha"] * x[j] / x[ideriv]^2 * grclgjk 
		res <- res - paramj["alpha"] / x[ideriv]^2 * clgjk * (j == mderiv) 
		res <- res + 2*paramj["alpha"] * x[j] / x[ideriv]^3 * clgjk * (ideriv == mderiv)
	}
	
	as.numeric((paramj["rmax"] - paramj["rmin"]) * res)
}




#______________
# premium difference based function

#examples of parameters c(rmin=0.02, rmax=0.98, mu=-2.890372, alpha=7.401976)
probj2kDiff <- function(x, j, k, paramj)
{
	
	n <- length(x)
	
	denom <- 1 + sum( exp(paramj["mu"] + paramj["alpha"] * (x[j] - x[-j])) )
	if(j == k)
	{
		num <- 1
		rtildemin <- 1 - paramj["rmax"] 
	}else
	{
		num <- exp(paramj["mu"] + paramj["alpha"] * (x[j] - x[k]) )
		rtildemin <- paramj["rmin"] / (n-1)
	}
	as.numeric( rtildemin + (paramj["rmax"] - paramj["rmin"]) * num / denom )
}



grprobj2kDiff <- function(x, j, k, ideriv, paramj)
{
	n <- length(x)
	clgjk <- probj2kDiff(x, j, k, paramj)

	if(ideriv == j)
	{
		res <- (-paramj["alpha"]) * sum( sapply( (1:n)[-j], function(l) probj2kDiff(x, j, l, paramj) ) ) * clgjk
	}else
	{
		res <- paramj["alpha"] * probj2kDiff(x, j, ideriv, paramj) * clgjk
	}
	
	if(j != k)
	{
		res <- res + (ideriv == j) * paramj["alpha"] * clgjk - (ideriv == k) * paramj["alpha"] * clgjk
	}
	
	as.numeric((paramj["rmax"] - paramj["rmin"]) * res)
}


grgrprobj2kDiff <- function(x, j, k, ideriv, mderiv, paramj)
{
	n <- length(x)
	clgji <- probj2kDiff(x, j, ideriv, paramj)
	clgjk <- probj2kDiff(x, j, k, paramj)
	grclgji <- grprobj2kDiff(x, j, ideriv, mderiv, paramj)
	grclgjk <- grprobj2kDiff(x, j, k, mderiv, paramj)
	
	if(ideriv == j)
	{
		sum1 <- sum( sapply( (1:n)[-j], function(l) grprobj2kDiff(x, j, l, mderiv, paramj) ) )
		sum2 <- sum( sapply( (1:n)[-j], function(l) probj2kDiff(x, j, l, paramj) ) )
		
		res <- -paramj["alpha"] * sum1 * clgjk - paramj["alpha"] * sum2 * grclgjk			
	}else
	{
		res <- paramj["alpha"] * grclgji * clgjk + paramj["alpha"] * clgji * grclgjk
	}	
	
	if(j != k && ideriv == j)
		res <- res + paramj["alpha"] * grclgjk 
	if(j != k && ideriv == k)
		res <- res - paramj["alpha"] * grclgjk 
	
	
	as.numeric((paramj["rmax"] - paramj["rmin"]) * res)
}



#__________________________________________
# TEST price ratio


Nplayer <- 3
n <- 10
xset <- matrix(runif(n*Nplayer, -5, 5), nrow=n, ncol=Nplayer)
xset <- rbind(xset, c(Inf, 1, 1))
xset <- rbind(xset, c(Inf, Inf, 1))
xset <- rbind(xset, c(-Inf, Inf, 1))

parlapse <- c(rmin=0.02, rmax=0.98, mu=-16.52114, alpha=13.63030)

res1 <- sapply(1:Nplayer, function(k)
	sapply(1:NROW(xset), function(l) 
		PROBj2kRatio(xset[l, ], 1, k, parlapse) )
	)


res2 <- sapply(1:Nplayer, function(k)
	sapply(1:NROW(xset), function(l) 
		probj2kRatio(xset[l, ], 1, k, parlapse) )
	)

res1 <- sapply(1:Nplayer, function(k)
	sapply(1:NROW(xset), function(l) 
		PROBj2kDiff(xset[l, ], 1, k, parlapse) )
	)


res2 <- sapply(1:Nplayer, function(k)
	sapply(1:NROW(xset), function(l) 
		probj2kDiff(xset[l, ], 1, k, parlapse) )
	)
