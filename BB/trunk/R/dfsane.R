  dfsane <- function(par, fn, M=5, maxiter=1000, lower=-Inf, upper=Inf, tol=1.e-07, trace=FALSE, ...) {
############################################################
# Non-monotone spectral method for finding a root of nonlinear systems
# LaCruz and Raydan M (2003): 
###############################################
# R translation (with minor modifications):  Ravi Varadhan, Johns Hopkins University, February 23, 2008.
################################################

######################################
#  local function
     lsm <- function(x, fn, F, fval, alfa, M, lastfv, eta, fcnt, bl, lower, upper, ...) {
	#  non-monotone line search of Grippo
######################################
	maxbl <- 100
	gamma <- 1.e-04
	sigma1 <- 0.1
	sigma2 <- 0.5
	lam1 <- 1
	lam2 <- 1
	cbl <- 0
	fmax <- max(lastfv)
	lsflag <- 1
	
#     Main loop
	while (cbl < maxbl & lsflag !=0) {

	d <- - alfa * F
	xnew <-  x + lam1 * d
	xnew[xnew < lower] <- lower[xnew < lower]
	xnew[xnew > upper] <- upper[xnew > upper]

	Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1

	if (class(Fnew) == "try-error" | any(is.nan(Fnew))) return(NULL) 
		else fune1 <- sum(Fnew * Fnew)

      if (fune1 <= (fmax + eta - (lam1^2*gamma*fval))) {
		if (cbl >= 1) bl <- bl + 1
          	lsflag <- 0
		fune <- fune1
		return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=lsflag, fune=fune))
	} 

	xnew <- x - lam2 * d
	xnew[xnew < lower] <- lower[xnew < lower]
	xnew[xnew > upper] <- upper[xnew > upper]

	Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1
	if (class(Fnew) == "try-error" | any(is.nan(Fnew))) return(NULL) 
		else fune2 <- sum(Fnew * Fnew)
      	if (fune2 <= (fmax + eta - (lam2^2*gamma*fval))) {
		if (cbl >= 1) bl <- bl + 1
          	lsflag <- 0
		fune <- fune2
		return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=lsflag, fune=fune))
	} 
#     Quadratic interpolation
      	lamc <- (2*fval*lam1^2) / (2 * (fune1 + (2*lam1 - 1)*fval))
		c1 <- sigma1 * lam1
		c2 <- sigma2 * lam1
		if (lamc < c1) lam1 <- c1
		else if (lamc > c2)lam1 <- c2
		else lam1 <- lamc

      	lamc <- (2*fval*lam2^2) / (2 * (fune2 + (2*lam2 - 1)*fval))
		c1 <- sigma1 * lam2
		c2 <- sigma2 * lam2
		if (lamc < c1) lam2 <- c1
		else if (lamc > c2)lam2 <- c2
		else lam2 <- lamc

		cbl <- cbl + 1
		lsflag <- 1
#		}
	}

	return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=lsflag, fune=fune))
	}

##########################################
#
#     Initialization
	n <- length(par)
	fcnt <- 0
	iter <- 0
	bl <- 0
	alfa <- 1
	eta <- 1
	ep <- 1.0e-10
	h <- 1.0e-07
	lastfv <- rep(0, M)
	stagn <- FALSE

  F <- try (fn(par, ...))

	if (class(F) == "try-error" | any(is.nan(F))) {
	cat(" Failure: Error in functional evaluation ")
	return(NULL)
	} else normF <- normFI <- sqrt(sum(F * F))

	if (trace) cat("||F(x0)||: ", normF, "\n")

      lastfv[1] <- normF^2
			
######################################
#     MAIN LOOP:  Iteration begins
######################################
	while (normF/sqrt(n) > tol & iter <= maxiter & !stagn) {

# Control of steplength
	if ((alfa <= ep) | (alfa >= 1/ep)) {
		if (normF > 1) alfa <- 1 
		else if (normF >= 1e-05 & normF <= 1)  alfa <- normF
		else if(normF < 1e-05)	alfa <- 1e-05
		}

#  non-monotone line search of Grippo
      ls.ret <-  lsm(x=par, fn=fn, F=F, fval=normF^2, alfa, M=M, lastfv=lastfv, eta, fcnt, bl, 
		lower, upper, ...)

	if(is.null(ls.ret)) {
	cat(" \n Failure: Error in functional evaluation \n")
	return(NULL)
	}

	flag <- ls.ret$lsflag

      if (flag == 1) {
	cat ("\n Failure: maximum number of steplength reductions is exceeded! \n")
	return(NULL)
	}

	Fnew <- ls.ret$Fnew
	pnew <- ls.ret$xnew
	fune <- ls.ret$fune
	fcnt <- ls.ret$fcnt
	bl <- ls.ret$bl

#     Calculate new steplength: alfa
	alfa <- sum((pnew - par)^2) / abs(sum((pnew - par) * (Fnew - F)))

	if (is.nan(alfa)) { 
		stagn <- TRUE
		next
		}
	par <- pnew
	F <- Fnew
	fun <- fune
	normF <- sqrt(fun)

	iter <- iter + 1
      	lastfv[ 1 + iter %% M] <- fun
	eta <- 1 / (2 ^iter)
	
      	if (trace) cat("\n iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")
	}
	if (normF/sqrt(n) <= tol) conv <- list(type=0, message="Successful convergence") 
	if (iter > maxiter) conv <- list(type=1, message="Maximum # iterations exceeded")
	if (stagn) conv <- list(type=2, message="Anomalous iteration")
	res <- normF / sqrt(length(par))

	return(list(par = par, resid = res, feval=fcnt, iter=iter, convergence=conv$type, message=conv$message))
      }

