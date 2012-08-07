
GNE.nseq.old <- function(xinit, Phi, jacPhi, argPhi=list(), argjac=list(), method="Newton", control=list(), ...)
{
	if(method == "default")
		method <- "Newton"
		
	nseq.old(xinit, Phi, jacPhi, argPhi, argjac, method, control, ...)	
}
