nseq.old <- function(xinit, Phi, jacPhi, argPhi, argjac, method, control, global="gline", 
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
		LM.param <- match.arg(con$LM.param, c("merit", "jacmerit", "min", "adaptive"))
		
		if(LM.param == "adaptive")
			test.try <- try( nseq.LM.adapt(xinit, wrapPhi, wrapJac, argPhi, argjac, 
							global=global, control=con), silent=silent)	

		if(LM.param != "adaptive")
			test.try <- try( nseq.LM(xinit, wrapPhi, wrapJac, argPhi=argPhi, argjac=argjac, 
							global=global, control=con), silent=silent)
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

