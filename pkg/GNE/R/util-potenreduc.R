#functions of the potential reduction algorithm



#logarithm potential function	
potential.ce <- function(u, n, zeta)
{
	if(any(u[-(1:n)] < 0))
		stop(paste("u has negative components:", paste(u, collapse=", "), "\n"))
	zeta * log( sum(u^2) ) - sum( log( u[-(1:n)] ) )
}

#gradient of logarithm potential function		
gradpotential.ce <- function(u, n, zeta)	
{
	normu <- sum(u^2)
	c( 2*zeta / normu * u[1:n], 2*zeta / normu * u[-(1:n)] - 1/u[-(1:n)] )
}	



#merit function for the constrained equation	
#z = (x, lam, w)
#with size (n, m, m)
psi.ce <- function(z, dimx, dimlam, Hfinal, argfun, zeta)
as.numeric( potential.ce( Hfinal(z, argfun=argfun), sum(dimx), zeta) )

#gradient of the merit function	
#z = (x, lam, w)
#with size (n, m, m)
gradpsi.ce <- function(z, dimx, dimlam, Hfinal, jacHfinal, 
argfun, argjac, zeta)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
	
	Hz <- Hfinal(z, argfun=argfun)
	Jz <- jacHfinal(z, argjac=argjac)
	
#t(Jz) =
#	t(Jz[1:n, 1:n])				t(Jz[(n+1):(n+m), 1:n])		0
#	t(Jz[1:n, (n+1):(n+m)])		0							diag(w)
#	0							Identity_m					diag(lam)			
	
	grpHz <- gradpotential.ce(Hz, n, zeta)
	
	c(	t(Jz[1:n, 1:n]) %*% grpHz[1:n] + t(Jz[(n+1):(n+m), 1:n]) %*% grpHz[(n+1):(n+m)], 
	  t(Jz[1:n, (n+1):(n+m)])  %*% grpHz[1:n] + diag(w) %*% grpHz[(n+m+1):(n+2*m)], 
	  grpHz[1:m+n] + diag(lam) %*% grpHz[(n+m+1):(n+2*m)]		)
}


#check interior
#z = (x, lam, w)
#with size (n, m, m)
checkint.ce <- function(z, dimx, dimlam)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
#	cat("checkint.ce\n")
#	print(cbind(lam, w))
	
	all(lam > 0, w > 0)
}

