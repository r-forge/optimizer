fixedpoint <- function(xinit, method=c("pure","UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"), 
	yfunc, stepfunc, Vfunc, control=list(), ...)
{
	method <- match.arg(method, c("pure","UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"))
	
	if(method == "UR")
		if(missing(yfunc) || missing(xinit) || missing(stepfunc) || missing(Vfunc))
			stop("missing parameters for UR method")
	
	if(method == "vH")
		if(missing(yfunc) || missing(xinit) || missing(Vfunc))
			stop("missing parameters for vH method")
	
	
	testmyfunc <- yfunc(xinit, ...)

	testVfunc <- Vfunc(xinit, ...)
	
	if(class(testmyfunc) %in% c("integer", "numeric"))
	{	
		inner.counts.fn <- inner.counts.vh <- NULL
		inner.iter.fn <- inner.iter.vh <- NULL
		noitercount <- TRUE
	}else if(class(testmyfunc) == "list")
	{
		noitercount <- FALSE
		
		inner.counts.fn <- c(0, 0)
		inner.iter.fn <- 0
		
		wrapyfunc <- function(x, ...) 
		{
			fx <- yfunc(x, ...)
			inner.counts.fn <<- inner.counts.fn + fx$counts
			inner.iter.fn <<- inner.iter.fn + fx$iter
			fx$par
		}
	
		if(class(testVfunc) == "list")
		{
			inner.counts.vh <- c(0, 0)
			inner.iter.vh <- 0
			wrapVfunc <- function(x, ...) 
			{
				vx <- Vfunc(x, ...)
				inner.counts.vh <<- inner.counts.vh + vx$counts
				inner.iter.vh <<- inner.iter.vh + vx$iter
				vx$value
			}			
		}else
		{
			inner.counts.vh <- NULL
			inner.iter.vh <- NULL			
		}	
		
	}else
		stop("Unsupported class for argument functions.")
		
	
	#default control parameters
	con <- list(sigma=9/10, beta=1/2, tol=1e-6, maxit=100, echo=TRUE)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	if(method == "pure")
	{
		if(noitercount)
			myGNE <- try( pureFP(xinit, yfunc, Vfunc, control=con, ...), silent=FALSE)
		else
			myGNE <- try( pureFP(xinit, wrapyfunc, wrapVfunc, control=con, ...), silent=FALSE)
		
		if(class(myGNE) == "try-error")
		{
			cat("\n")	
			stop("Pure fixed-point method does not converge.")
		}	
	}
	
	
	if(method == "UR")
	{
		if(noitercount)
			myGNE <- try( relaxationAlgoUR(xinit, stepfunc, yfunc, Vfunc, control=con, ...), silent=FALSE)
		else
			myGNE <- try( relaxationAlgoUR(xinit, stepfunc, wrapyfunc, wrapVfunc, control=con, ...), silent=FALSE)

		if(class(myGNE) == "try-error")
		{
			cat("\n")	
			stop("Relaxation algorithm UR for the fixed-point pb does not converge.")
		}	
	}
	
	if(method == "vH")
	{
		if(noitercount)
			myGNE <- try( relaxationAlgoVH(xinit, yfunc, Vfunc, control=con, ...), silent=FALSE)
		else
			myGNE <- try( relaxationAlgoVH(xinit, wrapyfunc, wrapVfunc, control=con, ...), silent=FALSE)

		if(class(myGNE) == "try-error")
		{
			cat("\n")	
			stop("Relaxation algorithm vH for the fixed-point pb does not converge.")
		}	
	}
	
	if(method == "RRE" || method == "MPE" || method == "SqRRE" || method == "SqMPE")
	{
		if(noitercount)
			myGNE <- try( extrapolFP(xinit, yfunc, Vfunc, control=con, method=method, ...), silent=FALSE)
		else
			myGNE <- try( extrapolFP(xinit, wrapyfunc, wrapVfunc, control=con, method=method, ...), silent=FALSE)
		
		if(class(myGNE) == "try-error")
		{
			cat("\n")	
			stop("Relaxation algorithm vH for the fixed-point pb does not converge.")
		}	
	}
	
#	extrapolFP <- function(xinit, yfunc, Vfunc, control, method, ...)
	
	
	k <- myGNE$iter
	xk <- myGNE$par
	f_xk <- myGNE$value
	counts <- myGNE$counts

	list(par=xk, outer.counts=counts, outer.iter=k, code=(k >= con$maxit)*1 + (abs(f_xk) > con$tol)*10, 
		 inner.iter.fn=inner.iter.fn, inner.iter.vh=inner.iter.vh, 
		 inner.counts.fn=inner.counts.fn, inner.counts.vh=inner.counts.vh)

	
}	


#pure fixed point iteration
pureFP <- function(xinit, yfunc, Vfunc, control, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	k <- 1
	xk_1 <- xinit
	xk <- yfunc(xk_1, ...)
	V_xk <- Vfunc( xk, ... )
	counts <- c(1, 1)
	names(counts) <- c("yfunc", "Vfunc")

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
			sqrt(sum( (xk - xk_1)^2 )),  "\n")
	
	
	while( (abs( V_xk ) > tol || sqrt( sum( (xk - xk_1)^2 ) ) > tol) && k < maxit) 
	{
		xk_1 <- xk
		k <- k+1
	
		xk <- yfunc(xk_1, ...)
		V_xk <- Vfunc( xk, ... )
		
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
				sqrt(sum( (xk - xk_1)^2 )),  "\n")
	
		counts <- counts+1
	}
	
	list(par = xk, value=V_xk , counts=counts, iter = k, code=(k >= maxit)*1 + (abs(V_xk) > tol)*10)
}


#Uryasev and Rubinstein
relaxationAlgoUR <- function(xinit, stepfunc, yfunc, Vfunc, control, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	k <- 1
	xk_1 <- xinit
	xk <- yfunc(xk_1, ...)
	V_xk <- Vfunc( xk, ... )
	
	counts <- c(1, 1)
	names(counts) <- c("yfunc", "stepfunc")
	
	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
			sqrt(sum( (xk - xk_1)^2 )),  "\n")

	while( (abs( V_xk ) > tol || sqrt( sum( (xk - xk_1)^2 ) ) > tol) && k < maxit)
	{
		xk_1 <- xk
		k <- k+1
		
		alphak <- stepfunc(k)
		xk <- (1-alphak) * xk_1 + alphak * yfunc(xk_1, ...)
		
		V_xk <- Vfunc( xk, ... )

		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
				sqrt(sum( (xk - xk_1)^2 )),  "\n")
		if(echo >= 3)
			cat("step size", alphak, "\n")
		
		counts <- counts+1
	}
	
	list(par = xk, value=V_xk , counts=counts, iter = k, 
		 code=(k >= maxit)*1 + (abs( V_xk ) > tol || sqrt( sum( (xk - xk_1)^2 ) ) > tol) * 10)
}


#von Heusinger 
relaxationAlgoVH <- function(xinit, yfunc, Vfunc, control, ...)
{
	sigma <- control$sigma
	beta <- control$beta
	tol <- control$tol
	echo <- control$echo
	maxit <- control$maxit
	
	k <- 0
	xk_1 <- xinit
	xk <- yfunc(xk_1, ...)
	V_xk <- Vfunc( xk, ... )
	
	counts <- c(1, 1)
	names(counts) <- c("yfunc", "Vfunc")

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||V(x_k)||", abs(V_xk),  "\n")

	
	if(echo) 
		cat("**** k", k, "\n", "||y(x_k) - x_k||", sqrt( sum( (xk - xk_1)^2 ) ),  "\n", "||V(x_k)||", abs(V_xk), "\n", "y(x_k)\t", xk, "\n")

	while( abs( V_xk ) > tol && k < maxit)
	{
		k <- k+1
		xk_1 <- xk
		
		dk <- yfunc(xk, ...) - xk
		normdk <- sqrt(sum(dk^2))
		V_xk <- Vfunc(xk, ...)
		
		
		for(l in 0:16)
		{
			tk <- beta^l
			
			if(echo >= 3)
			{
				cat(l, "\t", Vfunc(xk + tk * dk, ...) , "\t <= ")
				cat(V_xk - control$sigma * tk^2 * normdk^2, "?\t")
				cat("tk", tk, "\n")
			}
			
			# use remark 4.32 			
			if(Vfunc(xk + tk * dk, ...) <= V_xk - sigma * tk^2 * normdk^2)
				break
			
		}
				
		xk <- xk_1 + tk * dk
		
		V_xk <- Vfunc(xk, ...)

		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||V(x_k)||", abs(V_xk),  "\n")

		
		if(echo) 
			cat("\n\n**** k", k, "\n", "||y(x_k) - x_k||", sqrt( sum( (xk - xk_1)^2 ) ),  "\n", "||V(x_k)||", abs(V_xk), "\n", "y(x_k)\t", xk, "\n")		
		
		counts <- counts+1
		counts[2] <- counts[2]+l
	}	
	
	list(par = xk, value=V_xk, counts=counts, iter = k, code=(k >= maxit)*1 + (abs(V_xk) > tol)*10 )
}



#extrapolation method for fixed point iteration
extrapolFP <- function(xinit, yfunc, Vfunc, control, method, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	k <- 1
	xk_1 <- xinit
	xk <- yfunc(xk_1, ...)
	V_xk <- Vfunc(xk, ...)
	
	counts <- c(1, 1)
	names(counts) <- c("yfunc", "Vfunc")

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
			sqrt(sum( (xk - xk_1)^2 )),  "\n")
	
	
	while( (abs( V_xk ) > tol || sqrt( sum( (xk - xk_1)^2 ) ) > tol) && k < maxit) 
	{
		xk_1 <- xk
		k <- k+1
		
		xk <- yfunc(xk_1, ...)
		xkp1 <- yfunc(xk, ...)
		#RRE/MPE cycle of order 1, equivalent to Aitken Delta process when xk is univariate
		Delta1_xk <- xk - xk_1
		Delta2_xk <- xkp1 - 2*xk + xk_1

		if(method == "RRE" || method == "SqRRE")
			mystep <- crossprod(Delta2_xk, Delta1_xk) / crossprod(Delta2_xk, Delta2_xk)
		if(method == "MPE" || method == "SqMPE")
			mystep <- crossprod(Delta1_xk, Delta1_xk) / crossprod(Delta1_xk, Delta2_xk)

		if(is.nan(mystep) || is.infinite(mystep))
		{
			warning("Delta square xk is too small: Use the last iterate.")
			xk <- xkp1
		}else if(substr(method, 1, 2) == "Sq")
		{
			xk <- xk_1 - 2*Delta1_xk * mystep + Delta2_xk * mystep^2
		}else #simple version
		{
			xk <- xk_1 - Delta1_xk * mystep
		}
		
		V_xk <- Vfunc(xk, ...)
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||V(x_k)||", abs(V_xk),  "\n", "||y(x_k) - x_k||", 
				sqrt(sum( (xk - xk_1)^2 )),  "\n")
		if(echo >= 3)		
			cat(" Delta(xk)", Delta1_xk, "\n", "Delta^2(xk)", Delta2_xk, "\n")
		
		counts <- counts + c(2,1)
	}
	
	list(par = xk, value=V_xk , counts=counts, iter = k, 
		 code=(k >= maxit)*1 + (abs( V_xk ) > tol || sqrt( sum( (xk - xk_1)^2 ) ) > tol)*10 )
}

