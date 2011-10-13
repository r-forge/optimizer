nseq <- function(xinit, Phi, jacPhi, argPhi, argjac, method, control, global="gline", 
	silent=TRUE, ...)	
{
	
	
	method <- match.arg(method, c("Newton", "Broyden", "Levenberg-Marquardt"))
	global <- match.arg(global, c("dbldog", "pwldog", "qline", "gline", "none"))
	
	
	#default control parameters
	con <- list(ftol=1e-6, maxit=100, trace=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	#wrapped functions
	wrapPhi <- function(x, argPhi, argjac)
		evalwitharglist(Phi, x, argPhi) 
	
	wrapJac <- function(x, argPhi, argjac)
		evalwitharglist(jacPhi, x, argjac)
	
	#basic tests for Phi
	test.try <- try( wrapPhi(xinit, argPhi), silent=TRUE )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				message="Can't evalate Phi(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Phi(xinit) has infinite or NaN values.", fvec=NA) )

	#basic tests for JacPhi
	test.try <- try( wrapJac(xinit, argjac=argjac), silent=TRUE )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evaluate Jac Phi(xinit).", fvec=NA) )
	if(any(is.nan(test.try)) || any(is.infinite(test.try)) )
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Jac Phi(xinit) has infinite or NaN values.", fvec=NA) )
	
	if(method != "Levenberg-Marquardt")
		test.try <- try( nleqslv(xinit, wrapPhi, wrapJac, argPhi=argPhi, argjac=argjac,
			method = method, global = global, control=con, ...), silent=silent)
	
	if(method == "Levenberg-Marquardt")
	{
		test.try <- try( nseq.LM(xinit, wrapPhi, wrapJac, argPhi=argPhi, argjac=argjac,
			method="Levenberg-Marquardt", global=global, control=con), silent=silent)
	}	
	
	
	if(class(test.try) == "try-error")
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."), fvec=NA)
	if(class(test.try) != "try-error")
		res <- list(par = test.try$x, value = sqrt(sum( test.try$fvec^2 )), 
					counts = c(phicnt = test.try$nfcnt, jaccnt = test.try$njcnt), 
					iter = test.try$njcnt, code = test.try$termcd, 
					message = test.try$message, fvec= test.try$fvec)
	res
}	


nseq.LM <- function(xinit, Phi, jacPhi, argPhi, argjac, control, global="gline", silent=TRUE, ...)	
{
	global <- match.arg(global, c("gline", "none"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, btol=1e-2, maxit=100, trace=0, sigma=1/2, 
				echofile=NULL, echograph="NULL", delta=2, initlnsrch=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	xk_1 <- xinit
	xk <- xinit
	fk <- Phi(xinit, argPhi=argPhi)
	Jacfk <- jacPhi(xinit, argjac=argjac)
	iter <- 0
	inner.iter <- 0
	nfcnt <- 0
	njcnt <- 0
	termcd <- 0

	while(termcd == 0 && iter < con$maxit) 
	{
		lambdak <- sqrt( sum( fk^2 ) )^con$delta
		A <- crossprod( Jacfk, Jacfk ) + lambdak * diag( length(xk) )
		b <- -crossprod( Jacfk, fk )
		mycatch <- try( dk <- solve(A, b) , silent=silent)
		
		if(class(mycatch) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(mycatch), "singul")[[1]]) == 2 )
				termcd <- 6
			break
		}
		
		if(global == "none")
		{
			xkp1 <- xk + dk
			inner.iter <- 0
		}
		
		if(global == "gline")
		{
			normfk <- crossprod(fk)/2
			
			inner.iter <- con$initlnsrch
			
			slopek <- crossprod(Jacfk %*% dk, fk)
			
			minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )
			
			stepk <- 1
			
			while(stepk > minstep)
			{
				stepk <- con$sigma^inner.iter
				normfp <- crossprod(Phi(xk + stepk*dk, argPhi=argPhi))/2
				
				#traces in R console	
				if(con$trace >= 3)
				cat("i", inner.iter, "\tlambda", stepk, "\tright term\t", normfk + con$btol * stepk * slopek, "\tleft term", normfp, "\n")			
				
				cat("largest\t", max(abs(Phi(xk + stepk*dk, argPhi=argPhi))), "\n")	
				
				#check Armijo condition
				if(normfp <= normfk + con$btol * stepk * slopek)
				{
				#cat("Armijo satisfied\n")
					break
				}
				
				#cat("\tnorm stepk*dk\t", sqrt(sum(stepk*dk^2)), "\n")
				inner.iter <- inner.iter + 1	
				
			}
			if(stepk <= minstep)
			{	
				termcd <- 3
				break
			}else if(normfp <= normfk + con$btol * stepk * slopek) 
			{	
				xkp1 <- xk + stepk*dk
			}else
			{
				stop("internal error in nseq.LM function.")
			}
		}
		

		
		xk_1 <- xk
		xk <- xkp1
		
		fk_1 <- fk	
		fk <- Phi(xk, argPhi=argPhi)
		Jacfk <- jacPhi(xk, argjac=argjac)
		
		if(any(is.nan(Jacfk)) || any(is.nan(fk)))
			termcd <- -10
			
		
		iter <- iter+1
		nfcnt <- nfcnt + 1 + inner.iter - con$initlnsrch
		njcnt <- njcnt + 1
		
		#termination criterion, see Schnabel algo A7.2.3
		if(max( abs( fk ) ) <= con$ftol)
			termcd <- 1
		if(iter >= con$maxit)
			termcd <- 4
		if(max( abs(xk - xk_1) / abs(xk) ) <= con$xtol)
			termcd <- 2
		
	}
	
	message <- NA
	if(termcd == 1)
		message <- "Function criterion near zero"
	if(termcd == 2)
		message <- "x-values within tolerance `xtol'"
	if(termcd == 3)
		message <- "No better point found (algorithm has stalled)"
	if(termcd == 4)
		message <- "Iteration limit exceeded"
	if(termcd == 5)
		message <- "Jacobian is too ill-conditioned"
	if(termcd == 6)
		message <- "Jacobian is singular"
	if(termcd == -10)
		message <- "Analytical Jacobian most likely incorrect"
	
    	
	
	list(x= as.vector(xk), fvec=fk, nfcnt=nfcnt, njcnt=njcnt, iter=iter, 
		 termcd=termcd, message=message)
	
}
