sane <- function(par, fn, M=5, maxit=1000, lower=-Inf, upper=Inf, tol=1.e-07, trace=FALSE, ...) {
############################################################
# Non-monotone spectral method for finding a root of nonlinear systems
# LaCruz and Raydan M (2003): 
###############################################
# R translation (with minor modifications):  Ravi Varadhan, Johns Hopkins University, February 23, 2008.
################################################
#
######################################
#  local function
#  non-monotone line search of Grippo
     lineSearch <- function(x, fn, F, fval, dg, M, lastfv, sgn, lambda, 
                     fcnt, bl, lower, upper, ...) {
######################################
	maxbl <- 100
	gamma <- 1.e-04
	sigma1 <- 0.1
	sigma2 <- 0.5
	cbl <- 0
	fmax <- max(lastfv)
	gpd <- -2 * abs(dg)
	lsflag <- 1
	
#     Main loop
	while (cbl < maxbl & lsflag !=0) {

	xnew <-  x + lambda* sgn * F
	xnew[xnew < lower] <- lower[xnew < lower]
	xnew[xnew > upper] <- upper[xnew > upper]

	Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1

	if (class(Fnew) == "try-error" | any(is.nan(Fnew))) return(NULL) 
		else fune <- sum(Fnew * Fnew)

      if (fune <= (fmax + lambda*gpd*gamma)) {
		if (cbl >= 1) bl <- bl + 1
          	lsflag <- 0
	} else {
#     Quadratic interpolation
      	lamc <- -(gpd*lambda^2) / (2 * (fune - fval - lambda*gpd))
		c1 <- sigma1 * lambda
		c2 <- sigma2 * lambda
		if (lamc < c1) lambda <- c1
		else if (lamc > c2)lambda <- c2
		else lambda <- lamc
		cbl <- cbl + 1
		lsflag <- 1
		}
	}

	return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, lambda=lambda, bl=bl, lsflag=lsflag, fune=fune))
	}

##########################################

#     Initialization
	n <- length(par)
	fcnt <- 0
	iter <- 0
	bl <- 0
	alfa <- 1
	ep <- 1.0e-08
	h <- 1.0e-07
	lastfv <- rep(0, M)
	stagn <- FALSE

      F <- try (fn(par, ...))

	if (class(F) == "try-error" | any(is.nan(F))) {
	cat(" Failure: Error in functional evaluation ")
	return(NULL)
	} else normF <- sqrt(sum(F * F))

	if (trace) cat("||F(x0)||: ", normF, "\n")

      lastfv[1] <- normF^2
			
######################################
#     MAIN LOOP:  Iteration begins
######################################
	while (normF/sqrt(n) > tol & iter <= maxit & !stagn) {

# Calculate the gradient of the merit function ||F(X)||
	xa <- par + h*F

	Fa <- try (fn(xa, ...)) 
      fcnt <- fcnt + 1
		
	if (class(Fa) == "try-error" | any(is.nan(Fa))) {
		cat(" Failure: Error in functional evaluation \n")
		break
		} else 	dg <- (sum(F * Fa) - normF^2) / h

      if (abs(dg/normF^2) < ep) {
	stagn <- TRUE
	next
	}

# Control of steplength
	if ((alfa <= ep) | (alfa >= 1/ep)) {
		if (normF > 1) alfa <- 1 
		else if (normF >= 1e-05 & normF <= 1)  alfa <- normF
		else if(normF < 1e-05)	alfa <- 1e-05
		}

	if (dg > 0) sgn <- -1
	else sgn <- 1

	lambda <- 1/alfa

#  non-monotone line search of Grippo
      ls.ret <-  lineSearch(x=par, fn=fn, F=F, fval=normF^2, dg=dg, M=M, lastfv=lastfv, sgn, lambda, fcnt, bl, 
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

	lambda <- ls.ret$lambda
	Fnew <- ls.ret$Fnew
	pnew <- ls.ret$xnew
	fune <- ls.ret$fune
	fcnt <- ls.ret$fcnt
	bl <- ls.ret$bl

#     Calculate new steplength: alfa
	alfa <- sum((Fnew - F)^2) / (lambda * sum(F*(F - Fnew)))
	if (abs(alfa) > 1e05 | is.nan(alfa)) {
	alfa <- sum(F * Fnew)
	alfa <- (1 - (alfa/normF^2)) / lambda
	}
	
	par <- pnew
	F <- Fnew
	fun <- fune
	normF <- sqrt(fun)

	iter <- iter + 1
      	lastfv[ 1 + iter %% M] <- fun
	
      	if (trace) cat("\n iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")
	}
	if (normF/sqrt(n) <= tol) conv <- list(type=0, message="Successful convergence") 
	if (iter > maxit) conv <- list(type=1, message="Maximum # iterations exceeded")
	if (stagn) conv <- list(type=2, message="Anomalous iteration")
	res <- normF / sqrt(length(par))

	return(list(par = par, resid = res, feval=fcnt, iter=iter, convergence=conv$type, message=conv$message))
      }


