  dfsane <- function(par, fn, method=2, control=list(), ...) {
############################################################
# Non-monotone spectral method for solving large-scale nonlinear systems of equations
# Without gradient information
# W LaCruz, JM Martinez, and M Raydan (2006): 
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
#   Note that method = 2 is the "default" since it performed better than others in our numerical experiments on a variety of problems
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
     lsm <- function(x, fn, F, fval, alfa, M, lastfv, eta, fcnt, bl, ...) {
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
    
#     Main loop
    while (cbl < maxbl) {

    d <- - alfa * F
    xnew <-  x + lam1 * d

    Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1

    if (class(Fnew) == "try-error" | any(is.nan(Fnew)) | any(is.na(Fnew)) ) return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA))
        fune1 <- sum(Fnew * Fnew)

      if (fune1 <= (fmax + eta - (lam1^2*gamma*fval))) {
        if (cbl >= 1) bl <- bl + 1
        return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=0, fune=fune1))
    } 

    xnew <- x - lam2 * d

    Fnew <- try(fn(xnew, ...))
        fcnt = fcnt + 1
    if (class(Fnew) == "try-error" | any(is.nan(Fnew)) | any(is.na(Fnew)) ) return(list(xnew=NA, Fnew=NA, fcnt=fcnt, bl=bl, lsflag=1, fune=NA)) 
	  fune2 <- sum(Fnew * Fnew)

        if (fune2 <= (fmax + eta - (lam2^2*gamma*fval))) {
        if (cbl >= 1) bl <- bl + 1
        return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=0, fune=fune2))
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
    }

    return(list(xnew=xnew, Fnew=Fnew, fcnt=fcnt, bl=bl, lsflag=2, fune=fune))
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
    lastfv <- rep(0, M)
	
    kount <- 0
  
	F <- try (fn(par, ...))

    if (class(F) == "try-error" | any(is.nan(F))) {
    cat(" Failure: Error in initial functional evaluation \n Try another starting value \n")
    return(NULL)
    } else F0 <- normF <- sqrt(sum(F * F))

    if (trace) cat("||F(x0)||: ", F0, "\n")
	pbest <- par
	normF.best <- normF

      lastfv[1] <- normF^2
            
######################################
#     MAIN LOOP:  Iteration begins
######################################
    while (normF/sqrt(n) > tol & iter <= maxit) {

# Control of steplength
    if ((abs(alfa) <= ep) | (abs(alfa) >= 1/ep)) {
        if (normF > 1) alfa <- 1 
        else if (normF >= 1e-05 & normF <= 1)  alfa <- 1 / normF
        else if(normF < 1e-05)  alfa <- 1e05
        }

#  non-monotone line search of Grippo
      ls.ret <-  lsm(x=par, fn=fn, F=F, fval=normF^2, alfa, M=M, lastfv=lastfv, eta, fcnt, bl, ...)

    fcnt <- ls.ret$fcnt
    bl <- ls.ret$bl
    flag <- ls.ret$lsflag

      if (flag > 0) {
	par <- pbest
	normF <- normF.best
	break 
	}

    Fnew <- ls.ret$Fnew
    pnew <- ls.ret$xnew
    fune <- ls.ret$fune

#     Calculate new steplength: alfa
#   Three Barzilai-Borwein steplength options.  
#
    pF <- sum((pnew - par) * (Fnew - F))
    pp <- sum((pnew - par)^2)
    FF <- sum((Fnew - F)^2)
    if (method==1) alfa <- pp / pF
    if (method==2) alfa <- pF / FF 
    if (method==3) alfa <- sign(pF) * sqrt(pp / FF)

    if (is.nan(alfa)) alfa <- 1.e-10	

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
    eta <- F0 / (iter + 1)^2
    
    if (trace & (iter%%triter == 0)) cat("\n iteration: ",iter, " ||F(xn)|| =  ", normF, "\n")

    }   # End of main loop

    if (flag==0) {
    if (normF/sqrt(n) <= tol) conv <- list(type=0, message="Successful convergence") 
    if (iter > maxit) conv <- list(type=1, message="Maximum limit for iterations exceeded")
    }	
    if (flag==1) conv <- list(type=2, message="Failure: Error in function evaluation")
    if (flag==2) conv <- list(type=3, message="Failure: Maximum limit on steplength reductions exceeded")
    res <- normF / sqrt(length(par))
    
    absred <- abs(normF - F0)
    relred <- 1 - abs((normF - F0)/F0)

    return(list(par = par, resid = res, abs.reduction = absred, rel.reduction=relred, feval=fcnt, iter=iter, 
    convergence=conv$type, message=conv$message)) }
