NewtonKKT <- function(xinit, 
	method=c("Newton", "Levenberg-Marquardt","BFGS", "L-BFGS-B"), 
	phi, jacphi, control=list(), phi.arg=list(), jacphi.arg=list(), 
	extrapol=c("none", "RRE", "MPE", "SqRRE", "SqMPE"), bounds,
	silent=TRUE)
{
	cat("Depreciated function.\n")
	
	method <- match.arg(method, c("Newton", "Levenberg-Marquardt","BFGS", "L-BFGS-B"))
	extrapol <- match.arg(extrapol, c("none", "RRE", "MPE", "SqRRE", "SqMPE"))
	
	#default control parameters
	con <- list(tol=1e-6, maxit=100, echo=0)
	namc <- names(con)
	con[namc <- names(control)] <- control

	if(method %in% c("Newton", "Levenberg-Marquardt"))
		myGNE <- try( EKKT_by_NLM(xinit, phi, jacphi, control=con, phi.arg, jacphi.arg, 
								  method, extrapol), silent=silent)
	if(method %in% c("BFGS", "L-BFGS-B"))
		myGNE <- try( EKKT_by_BFGS(xinit, phi, jacphi, control=con, phi.arg, jacphi.arg, 
								   method, extrapol, bounds), silent=silent)
	
	if(class(myGNE) == "try-error")
	{
		cat("\n")
		stop("Newton type method to solve the KKT system does not converge.")
	}else
		myGNE
}



#KKT system solving method by Newton or Levenberg-Marquardt method with possibly MPE/RRE cycle
EKKT_by_NLM <- function(xinit, phi, jacphi, control, phi.arg, jacphi.arg, method, extrapol)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	nbparam <- length(xinit)
	
	xk <- xinit
	xk_1 <- rep(Inf, nbparam)
	
	phik <- evalwitharglist(phi, xk, phi.arg)
	Jphik <- evalwitharglist(jacphi, xk, jacphi.arg)
	
	meritfunc <- function(x)
	{
		res <- sqrt( sum( evalwitharglist(phi, x, phi.arg)^2 ) )
#		cat("\nres merit", res, "\n")
		res
	}
	
	
	iter <- 0
	counts <- c(1,1)
	names(counts) <- c("phi", "jacphi")

	if(echo >= 1)
		cat("**** k", iter, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||phi(x_k)||", sqrt( sum(phik^2) ),  "\n")
	if(echo >= 3)
	{
		cat(" Jac Phi(x_k)\n")
		print(Jphik)
	}
	
	
	if(extrapol == "none")
	while(sqrt( sum( phik^2 ) ) > tol && iter < maxit) #sqrt( sum( (xk-xk_1)^2 ) ) > tol
	{
		#next step				
		if(method == "Newton")
			xkp1 <- NewtonNext(xk, phik, Jphik)
		if(method == "Levenberg-Marquardt")
			xkp1 <- LevenMarqNext(xk, phik, Jphik)
		
		#update
		xk_1 <- xk
		xk <- xkp1
		
		phik <- evalwitharglist(phi, xk, phi.arg)
		Jphik <- evalwitharglist(jacphi, xk, jacphi.arg)
		
		iter <- iter+1

		if(echo >= 1)
			cat("**** k", iter, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||phi(x_k)||", sqrt( sum(phik^2) ),  "\n")
		if(echo >= 3)
		{
			cat(" Jac Phi(x_k)\n")
			print(Jphik)
		}

	}

	if(extrapol %in% c("RRE", "MPE", "SqRRE", "SqMPE"))
	while(sqrt( sum( phik^2 ) ) > tol && iter < maxit) 
	{
		#next step				
		if(method == "Newton")
			xkp1 <- NewtonNext(xk, phik, Jphik)
		if(method == "Levenberg-Marquardt")
			xkp1 <- LevenMarqNext(xk, phik, Jphik)
		
		phikp1 <- evalwitharglist(phi, xkp1, phi.arg)
		Jphikp1 <- evalwitharglist(jacphi, xkp1, jacphi.arg)
		
		if(method == "Newton")
			xkp2 <- NewtonNext(xkp1, phikp1, Jphikp1)
		if(method == "Levenberg-Marquardt")
			xkp2 <- LevenMarqNext(xkp1, phikp1, Jphikp1)
		
		#RRE/MPE cycle of order 1
		Delta1_xk <- xkp1 - xk
		Delta2_xk <- xkp2 - 2*xkp1 + xk
		if(extrapol == "RRE" || extrapol == "SqRRE")
			mystep <- crossprod(Delta2_xk, Delta1_xk) / crossprod(Delta2_xk, Delta2_xk)
		if(extrapol == "MPE" || extrapol == "SqMPE")
			mystep <- crossprod(Delta1_xk, Delta1_xk) / crossprod(Delta1_xk, Delta2_xk)


		
		if(is.nan(mystep) || is.infinite(mystep))
		{
			warning("Delta square xk is too small: Use the last iterate.")
			tk <- xkp2
		}else if(substr(extrapol, 1, 2) == "Sq")
		{
			tk <- xk_1 - 2*Delta1_xk * mystep + Delta2_xk * mystep^2
		}else #simple version
		{
			tk <- xk - (xkp1 - xk) * mystep
		}
		
		#update
		xk_1 <- xk
		xk <- tk
		
		phik <- evalwitharglist(phi, xk, phi.arg)
		Jphik <- evalwitharglist(jacphi, xk, jacphi.arg)
		
		iter <- iter+1
		
		if(echo >= 1)
			cat("**** k", iter, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||phi(x_k)||", sqrt( sum(phik^2) ),  "\n")
		if(echo >= 3)
		{
			cat(" Jac Phi(x_k)\n")
			print(Jphik)
		}
		
	}	
	
	if(extrapol == "none")
		counts <- counts + iter
	if(extrapol %in% c("Aitken", "RRE", "MPE"))
		counts <- counts + 2*iter
		
	list(par= as.vector(xk), value=sqrt( sum( phik^2 ) ), counts=counts, iter=iter, 
		 code=(iter >= maxit)*1, message=NULL)
}


#KKT system solving method by BFGS/L-BFGS-B method
EKKT_by_BFGS <- function(xinit, phi, jacphi, control, phi.arg, jacphi.arg, method, extrapol, bounds)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	nbparam <- length(xinit)
	
	xk <- xinit
	
	meritfunc <- function(x)
	{
		res <- sqrt( sum( evalwitharglist(phi, x, phi.arg)^2 ) )
#		cat("\nres merit", res, "\n")
		res
	}
	
	gradmeritfunc <- function(x)
	{
		phix <- evalwitharglist(phi, x, phi.arg)	
		res <- crossprod( evalwitharglist(jacphi, x, jacphi.arg) , phix ) / sqrt( sum( phix^2 ) )
#		cat("grad merit", res, "\n\n")
		res
	}
	
	testinit <- meritfunc(xk)

	
	testinit <- gradmeritfunc(xk)

	
	if(method == "BFGS")
		resoptim <- optim(xinit, meritfunc, gradmeritfunc, method="BFGS", 
					  control=list(maxit=maxit, trace=echo, REPORT=1, abstol=tol))
	if(method == "L-BFGS-B")
		resoptim <- optim(xinit, meritfunc, gradmeritfunc, method="L-BFGS-B", lower=bounds[1], upper=bounds[2],
					  control=list(maxit=maxit, trace=echo, REPORT=1, pgtol=tol))
	
	if(echo >= 2)
	{
		print(resoptim)
		cat("\n")
	}
	
	k <- as.numeric(resoptim$counts[2])-1
	xk <- as.numeric(resoptim$par)
	f_xk <- as.numeric(resoptim$value)
	
	iter <- k
	counts <- resoptim$counts
	counts[1] <- sum(counts)
	names(counts) <- c("phi", "jacphi")
	
	
	list(par=xk, value=f_xk, counts=counts, iter=iter, code=(iter >= maxit)*1 + (abs(f_xk) > tol)*1,
		message=resoptim$message)
}


