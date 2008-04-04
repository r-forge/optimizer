sane <- function(par, fn, method=1, control=list(), ...) {
############################################################
# Non-monotone spectral method for finding a root of nonlinear systems
# LaCruz and Raydan M (Optimization Methods and Software 2003): 
###############################################
#
# R adaptation, with significant modifications, by Ravi Varadhan, Johns Hopkins University, March 25, 2008.
#
#   Most important modification is the availability of different options for Barzilai-Borwein steplengths
#   Three different Barzilai-Borwein steplength options can be chosen.
#   method = 1 is the steplength used in LaCruz, Martinez, and Raydan (2006)  
#   method = 2 is another BB steplength proposed in Barzilai and Borwein's (1988) original paper 
#   method = 3, is a new steplength, first proposed in Varadhan and Roland (2008).
#
#   Note that method = 1 is the "default" ; however, there doesn't seem to be much difference between 3 schemes
#
# Please refer to Varadhan and Gilbert (2008, unpublished) for details
#
################################################
    ctrl <- list(M=10, maxit=1500, tol=1e-07, trace=TRUE, triter=10) # defaults
    ctrl[names(control)] <- control
    M     <- ctrl$M
    maxit <- ctrl$maxit
    tol   <- ctrl$tol
    trace <- ctrl$trace
    triter <- ctrl$triter

######################################
#  local function
#  non-monotone line search of Grippo
     lineSearch <- function(x, fn, F, fval, dg, M, lastfv, sgn, lambda, 
                     fcnt, bl, ...) {
######################################
    maxbl <- 100
    gamma <- 1.e-04
    sigma1 <- 0.1
    sigma2 <- 0.5
    cbl <- 0
    fmax <- max(lastfv)
    gpd <- -2 * abs(dg)
    
#     Main loop
    while (cbl < maxbl) {

    xnew <-  x + lambda* sgn * F

    Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1

    if (class(Fnew) == "try-error" | any(is.nan(Fnew)) | any(is.na(Fnew)) | any(is.infinite(Fnew))) return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA))
        else fune <- sum(Fnew * Fnew)

      if (fune <= (fmax + lambda*gpd*gamma)) {
        if (cbl >= 1) bl <- bl + 1
	return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, lambda=lambda, bl=bl, lsflag=0, fune=fune))
    } else {
#     Quadratic interpolation
        lamc <- -(gpd*lambda^2) / (2 * (fune - fval - lambda*gpd))
        c1 <- sigma1 * lambda
        c2 <- sigma2 * lambda
        if (lamc < c1) lambda <- c1
        else if (lamc > c2)lambda <- c2
        else lambda <- lamc
        cbl <- cbl + 1
        }
    }
    return(list(xnew=NA, Fnew=NA, fcnt=fcnt, lambda=NA, bl=bl, lsflag=2, fune=NA))
    }

##########################################

#     Initialization
    n <- length(par)
    fcnt <- 0
    iter <- 0
    bl <- 0
    alfa <- 1
    eps <- 1.0e-08
    h <- 1.0e-07
    lastfv <- rep(0, M)
    stagn <- FALSE

      F <- try (fn(par, ...))

    if (class(F) == "try-error" | any(is.nan(F)) | any(is.na(F)) ) {
    cat(" Failure: Error in initial functional evaluation ")
    return(NULL)
    } else F0 <- normF <- sqrt(sum(F * F))

    if (trace) cat("||F(x0)||: ", normF, "\n")

	pbest <- par
	normF.best <- normF
      lastfv[1] <- normF^2
            
######################################
#     MAIN LOOP:  Iteration begins
######################################
    while (normF/sqrt(n) > tol & iter <= maxit & !stagn) {

# Calculate the gradient of the merit function ||F(X)||
    xa <- par + h*F

    Fa <- try (fn(xa, ...)) 
      fcnt <- fcnt + 1
        
    if (class(Fa) == "try-error" | any(is.nan(Fa)) | any(is.na(Fa)) ) {
		flag <- 1
		break
        } 

	dg <- (sum(F * Fa) - normF^2) / h
     if (abs(dg/normF^2) < eps) {
	flag <- 3
	break
	}

# Control of steplength
    if ((alfa <= eps) | (alfa >= 1/eps)) {
        if (normF > 1) alfa <- 1 
        else if (normF >= 1e-05 & normF <= 1)  alfa <- normF
        else if(normF < 1e-05)  alfa <- 1e-05
        }

    if (dg > 0) sgn <- -1
    else sgn <- 1

    lambda <- 1/alfa

#  non-monotone line search of Grippo
      ls.ret <-  lineSearch(x=par, fn=fn, F=F, fval=normF^2, dg=dg, M=M, lastfv=lastfv, sgn, lambda, fcnt, bl, ...)

    fcnt <- ls.ret$fcnt
    bl <- ls.ret$bl
    flag <- ls.ret$lsflag

      if (flag > 0) break 

    lambda <- ls.ret$lambda
    Fnew <- ls.ret$Fnew
    pnew <- ls.ret$xnew
    fune <- ls.ret$fune

#     Calculate new steplength: alfa

    if (method==1) alfa <- sum(F*(F - Fnew)) / (lambda * sum(F*F))
    if (method==2) alfa <- sum((F - Fnew)^2) / (lambda * sum(F*(F - Fnew)))
    if (method==3) alfa <- sign(sum(F*(F - Fnew))) * sqrt(sum((F - Fnew)^2) / (lambda^2 * sum(F*F)))

    if (is.nan(alfa)) alfa <- 1.e-08
    
    par <- pnew
    F <- Fnew
    fun <- fune
    normF <- sqrt(fun)

	if (normF < normF.best) {
	pbest <- par
	normF.best <- normF
	}

    iter <- iter + 1
        lastfv[ 1 + iter %% M] <- fun
    
    if (trace & (iter%%triter == 0)) cat("\n iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")

    }  # End of main loop

	if (flag==0) {
    if (normF/sqrt(n) <= tol) conv <- list(type=0, message="Successful convergence") 
    if (iter > maxit) conv <- list(type=1, message="Maximum number of iterations exceeded")
	} else {
	par <- pbest
	normF <- normF.best
    if (flag==1) conv <- list(type=2, message="Failure: Error in function evaluation")
    if (flag==2) conv <- list(type=3, message="Failure: Maximum limit on steplength reductions exceeded")
    if (flag==3) conv <- list(type=4, message="Failure: Anomalous iteration")
	}
    res <- normF / sqrt(length(par))

    absred <- abs(normF - F0)
    relred <- 1 - abs((normF - F0)/F0)

    return(list(par = par, resid = res, abs.reduction = absred, rel.reduction=relred, feval=fcnt, iter=iter, 
    convergence=conv$type, message=conv$message)) 

      }


