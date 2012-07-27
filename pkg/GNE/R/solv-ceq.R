ceq <- function(xinit, dimx, dimlam, dimw, 
	Hfinal, jacHfinal, argH, argjac,  
	method, control, global="gline", silent=TRUE, ...)	
{
	method <- match.arg(method, "IP")
	global <- match.arg(global, c("qline", "gline", "none"))


	#default control parameters
	con <- list(ftol=1e-6, maxit=100, trace=0)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	test.try <- try( ceq.IP(xinit, dimx, dimlam, dimw, Hfinal, jacHfinal, 
		argH, argjac, merit=psi.ce, gradmerit=gradpsi.ce, checkint=checkint.ce, 
		control= control, global=global, silent=silent), silent=silent)		
	
	if(class(test.try) == "try-error")
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."), 
					fvec=NA)
	if(class(test.try) != "try-error")
		res <- list(par = test.try$x, value = sqrt(sum( test.try$fvec^2 )), 
					counts = c(phicnt = test.try$nfcnt, jaccnt = test.try$njcnt), 
					iter = test.try$njcnt, code = test.try$termcd, 
					message = test.try$message, fvec= test.try$fvec)
	res
}	


ceq.IP <- function(xinit, dimx, dimlam, dimw, Hfinal, jacHfinal, argH, argjac,
	merit, gradmerit, checkint, control, global, silent=TRUE)	
{	
	global <- match.arg(global, c("gline", "qline", "none"))
	
	#default control parameters
	con <- list(ftol=1e-6, xtol=1e-5, btol=1e-2, maxit=100, trace=0, sigma=1/2, 
				echofile=NULL, delta=1, forcingpar=.1, zeta=length(xinit)/2)
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	xk_1 <- xinit
	xk <- xinit
	fk <- Hfinal(xinit, argH=argH)
	sigmak <- con$forcingpar
	
	iter <- 0
	inner.iter <- 0
	nfcnt <- 0
	njcnt <- 0
	termcd <- 0
	
	#traces in R console
	if(con$trace >= 1 && is.null(con$echofile))
		cat("**** k", iter, "\n")
	if(con$trace >= 1 && is.null(con$echofile))
		cat("x_k", xk, "\n")	
	if(con$trace >= 2 && is.null(con$echofile))
		cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
	

	while(termcd == 0 && iter < con$maxit) 
	{
		
		dk <- ceq.IP.direction(xk, dimx, dimlam, dimw, Hfinal, jacHfinal, 
			argH, argjac, sigmak, silent=silent)
		
		if(class(dk) == "try-error")
		{
			termcd <- 5
			if( length(strsplit(as.character(dk), "singul")[[1]]) == 2 )
				termcd <- 6
			break
		}
		
		if(global == "none")
		{
			xkp1 <- xk + dk
			inner.iter <- 0
		}
		
		if(global != "none")
		{	
			inner.iter <- 0
			minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )
			
			slopek <- crossprod( gradmerit(xk, dimx, dimlam, dimw, 
				Hfinal, jacHfinal, argH, argjac, con$zeta), dk)
			stepk <- 1			
			
			if(global == "gline")
				LSres <- try( linesearch.geom.cond(xk, dk, slopek, con, merit= merit, 
					checkint=checkint, dimx=dimx, dimlam=dimlam, dimw=dimw, 
					Hfinal=Hfinal, argH=argH, zeta=con$zeta) )
			
			
			if(global == "qline")
				LSres <- try( linesearch.quad.cond(xk, dk, slopek, con, merit= merit, 
					checkint=checkint, dimx=dimx, dimlam=dimlam, dimw=dimw, 
					Hfinal=Hfinal, argH=argH, zeta=con$zeta) )
			if(class(LSres) == "try-error")
				stop("internal error in line search function.")
			
			if(LSres$stepk <= minstep)
			{	
				termcd <- 3
				break
			}else if(LSres$normfp <= LSres$normfk + con$btol * LSres$stepk * slopek) 
			{	
				xkp1 <- xk + LSres$stepk*dk
			}else
			{
				stop("internal error in ceq.IP function.")
			}
		}
		

		
		xk_1 <- xk
		xk <- xkp1
		
		fk_1 <- fk	
		fk <- Hfinal(xk, argH=argH)		

		if(any(is.nan(fk)))
			termcd <- -10
			
		
		iter <- iter+1
		nfcnt <- nfcnt + 1 + LSres$inner.iter 
		njcnt <- njcnt + 1

		#traces in R console
		if(con$trace >= 1 && is.null(con$echofile))
			cat("**** k", iter, "\n")		
		if(con$trace >= 1 && is.null(con$echofile))
			cat("x_k", xk, "\n")
		if(con$trace >= 2 && is.null(con$echofile))
			cat(" ||f(x_k)||", sqrt( sum(fk^2) ),  "\n")
		
		
		
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


#z = (x, lam, w)
#with size (n, m, m)
ceq.IP.direction <- function(z, dimx, dimlam, dimw, Hfinal, jacHfinal, 
	argH, argjac, sigma, silent=TRUE)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
		
	Hz <- Hfinal(z, argH=argH)
	Jz <- jacHfinal(z, argjac=argjac)
	
	a <- c( rep(0, n), rep(1, 2*m) )
	
	btot <- sigma * sum(a*Hz) / sum(a^2) * a - Hz
	bx <- btot[1:n]
	blam <- btot[(n+1):(n+m)]
	bw <- btot[(n+m+1):(n+2*m)]
	
	Ex <- Jz[1:n, (n+1):(n+m)] 
	Jacgx <- Jz[(n+1):(n+m), 1:n]
	
	mat4x <- Jz[1:n, 1:n] + Ex %*% diag(lam/w) %*% Jacgx
	vec4x <- bx - Ex %*% diag(1/w) %*% bw + Ex %*% diag(lam/w) %*% blam
	
	d4x <- try( qr.solve(mat4x, vec4x), silent=silent)
	if(class(d4x) != "try-error")
	{
		d4w <- blam - Jacgx %*% d4x
		d4lam <- diag(1/w) %*% (bw - diag(lam) %*% d4w)
		return(c(d4x, d4lam, d4w))
	}else
		return(d4x)	
}




#logarithm potential function	
potential.ce <- function(u, n, zeta)
{
	if(any(u[-(1:n)] < 0))
		stop(paste("u has negative components:", u, "\n"))
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
psi.ce <- function(z, dimx, dimlam, dimw, Hfinal, argH, zeta)
	as.numeric( potential.ce( Hfinal(z, argH=argH), sum(dimx), zeta) )

#gradient of the merit function	
#z = (x, lam, w)
#with size (n, m, m)
gradpsi.ce <- function(z, dimx, dimlam, dimw, Hfinal, jacHfinal, 
	argH, argjac, zeta)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
	
	Hz <- Hfinal(z, argH=argH)
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
checkint.ce <- function(z, dimx, dimlam, dimw)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]

	all(lam > 0, w > 0)
}

