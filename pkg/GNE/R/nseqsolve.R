nseq <- function(xinit, Phi, jacPhi, argPhi, argjac, method, control, global="gline", 
	silent=TRUE, ...)	
{
	
	method <- match.arg(method, c("Newton", "Broyden", "Levenberg-Marquardt"))
	global <- match.arg(global, c("dbldog", "pwldog", "qline", "gline", "none"))

	#default control parameters
	con <- list(ftol=1e-6, maxit=100, trace=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	wrapPhi <- function(x, argphi, argjac)
		evalwitharglist(Phi, x, argPhi) 
	
	wrapJac <- function(x, argphi, argjac)
		evalwitharglist(jacPhi, x, argjac)
	

	test.try <- try( nleqslv(xinit, wrapPhi, wrapJac, argphi=argphi, argjac=argjac,
		method = method, global = global, control=con, ...), silent=silent)
	
	
	
	if(class(test.try) == "try-error")
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."))
	if(class(test.try) != "try-error")
		res <- list(par = test.try$x, value = sqrt(sum( test.try$fvec^2 )), 
					counts = c(phicnt = test.try$nfcnt, jaccnt = test.try$njcnt), 
					iter = test.try$njcnt, code = test.try$termcd, 
					message = test.try$message)
	res
}	

