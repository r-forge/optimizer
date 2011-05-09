
projector <- function(z, g, jacg, bounds=c(0, 10), echo=FALSE, ...)
{
	merit <- function(y, z) 1/2 * sum( (y - z)^2 )
	grmerit <- function(y, z) y - z	
	
	zinit <- rejection(g, length(z), bounds[1], bounds[2], ...)
	
	
	if(echo)
		cat("init", zinit,"\nz", z, "\n")
	

#	print(merit(zinit, z))
#	print(grmerit(zinit, z))
#	print(-g(zinit, ...))
#	print(-jacg(zinit, ...))
	
	res <- constrOptim.nl(zinit, 
				fn=function(y, ...) merit(y, z), 
				gr=function(y, ...) grmerit(y, z), 
				hin=function(y, ...) -g(y, ...), 
				hin.jac=function(y, ...) -jacg(y, ...), 
				control.outer=list(trace=FALSE),
				control.optim=list(trace=FALSE), ...)
	
	if(echo)
		print(res)
	
	res$par
}

