#H for solving the constrained equation
#z = (x, lam, w)
#with size (n, m, m)
Htemplate <- function(z, argH)
{
	x <- z[argH$dimz[1,"x"]:argH$dimz[2,"x"]]
	lam <- z[argH$dimz[1,"lam"]:argH$dimz[2,"lam"]]
	w <- z[argH$dimz[1,"w"]:argH$dimz[2,"w"]]
	
	c(	argH$lagreq(x, lam),
		argH$constr(x) + w,
		lam * w )
}

#z = (x, lam, w)
#with size (n, m, m)
jacHtemplate <- function(z, argjacH)
{
	x <- z[argjacH$dimz[1,"x"]:argjacH$dimz[2,"x"]]
	lam <- z[argjacH$dimz[1,"lam"]:argjacH$dimz[2,"lam"]]
	w <- z[argjacH$dimz[1,"w"]:argjacH$dimz[2,"w"]]
	
	n <- length(x)
	m <- length(w)
	
	rbind(
		cbind(argjacH$jaclagreq(x, lam), argjacH$diaggrconstr(x), diag(m)*0),
		cbind(argjacH$jacconstr(x), diag(m)*0, diag(m)),
		cbind(diag(m)*0, diag(w), diag(lam))				
	)
}


