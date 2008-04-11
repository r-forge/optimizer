spg <- function(par, fn, gr=NULL, method=3, project=NULL, 
           lower=-Inf, upper=Inf,  control=list(),  ... ) {

  # control defaults
  ctrl <- list(M=10, maxit=1500, gtol=1.e-05, maxfeval=10000, maximize=FALSE, 
        trace=TRUE, triter=10, grad.method="simple", eps=1.e-07) 
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])     

  ctrl[namc ] <- control
  M	<- ctrl$M
  maxit <- ctrl$maxit
  gtol  <- ctrl$gtol
  maxfeval <- ctrl$maxfeval
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  triter <- ctrl$triter
  grad.method <- ctrl$grad.method
  eps <- ctrl$eps
  ################ local function
  nmls <- function(p, f, d, gtd, lastfv, feval, func, maxfeval, ... ){
    # Non-monotone line search of Grippo with safe-guarded quadratic interpolation
    gamma <- 1.e-04
    fmax <- max(lastfv)
    alpha <- 1
    pnew <- p + alpha*d
    fnew <- try(func(pnew , ...),silent=TRUE)
    feval <- feval + 1
 
    if (class(fnew)=="try-error" | is.nan(fnew) )
            return(list(p=NA, f=NA, feval=NA, lsflag=1))
 
    while(fnew > fmax + gamma*alpha*gtd) {
    	if (alpha <= 0.1) alpha <- alpha/2
    	else {
    	    atemp <- -(gtd*alpha^2) / (2*(fnew - f - alpha*gtd))
    	    if (atemp < 0.1 | atemp > 0.9*alpha) atemp <- alpha/2
    	    alpha <- atemp
    	    }

    	pnew <- p + alpha*d
    	fnew <- try(func(pnew, ... ), silent=TRUE)
    	feval <- feval + 1
 
    	if (class(fnew)=="try-error" | is.nan(fnew) )
	       return(list(p=NA, f=NA, feval=NA, lsflag=1))
    	if (feval > maxfeval)
	       return(list(p=NA, f=NA, feval=NA, lsflag=2))
 
    	}  #while condition loop ends
 
    return(list(p=pnew, f=fnew, feval=feval, lsflag=0))
    }
 
  #############################################
  if (is.null(project)) project.box <- function(x, lower, upper) {
       # local function
       # Defined only when user-doesn't specify his/her own projection algorithm
       # Projecting to ensure that box-constraints are satisfied
       x[x < lower] <- lower[x < lower]
       x[x > upper] <- upper[x > upper]
       return(x)
       }
  #############################################

  #  Initialization
  lmin <- 1.e-30
  lmax <- 1.e30
  iter <-  feval <-  geval <- 0
  lastfv <- rep(-1.e99, M)
  fbest <- NA
 
  func <- if (maximize) function(x, ...) {-1 * fn(x, ...)}
                   else function(x, ...) fn(x, ...)

  # Project initial guess
  par <- if (is.null(project)) try(project.box(par, lower, upper), silent=TRUE)
                   else        try(project(par, ...), silent=TRUE)
 
  if (class(par) == "try-error" | any(is.nan(par)) | any(is.na(par)) ) {
     stop("Failure in projecting initial guess!")
  } else pbest <- par
 
  f <- try(func(par, ...),silent=TRUE)
  feval <- feval + 1

  if (class(f)=="try-error" | is.nan(f) | is.infinite(f) | is.na(f) )
    stop("Failure in initial function evaluation!")
    
  f0 <- fbest <- f
 
  if (is.null(gr)) require(numDeriv)
 
  g <- if (is.null(gr)) try(grad(func=func, x=par, method=grad.method, method.args=list(eps=eps, r=2), ...),silent=TRUE)
             else       try(gr(par, ...),silent=TRUE)
  
  geval <- geval + 1
 
  if (class(g)=="try-error" | any(is.nan(g))  ) 
    stop("Failure in initial gradient evaluation!")
 
  lastfv[1] <- f
  fbest     <- f
  pg   <- par - g
 
  pg <- if (is.null(project)) project.box(pg, lower, upper)
                   else       project(pg, ...)
 
  if (class(pg)=="try-error" | any(is.nan(pg))  ) 
    stop("Failure in initial projection!")
 
  pg <- pg - par

  pg2n <- sqrt(sum(pg*pg))
  pginfn <- max(abs(pg))
  gbest <- pg2n
  if (pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))
 
  if (trace) cat("iter: ",0, " f-value: ", f0, " pgrad: ",pginfn, "\n")

  #######################
  #  Main iterative loop
  #######################
  while( pginfn > gtol & iter <= maxit ) {
      iter <- iter + 1
      d <- par - lambda * g
 
      d <- if (is.null(project)) try(project.box(d, lower, upper), silent=TRUE)
                          else   try(project(d, ...), silent=TRUE)
 
      if (class(d) == "try-error" | any(is.nan(d))  ) {
        lsflag <- 4
        break
        }
 
      d <- d - par
      gtd <- sum(g * d)
 
      if(is.infinite(gtd)){
        lsflag <- 4
        break
        }
 
      nmls.ans <- nmls(par, f, d, gtd, lastfv, feval , func, maxfeval, ...)
      lsflag <- nmls.ans$lsflag
 
      if(lsflag != 0) break
 
      f <- nmls.ans$f
      pnew <- nmls.ans$p
      feval <- nmls.ans$feval
      lastfv[(iter %% M) + 1] <- f
 
      gnew <- if (is.null(gr)) try(grad(func=func, x=pnew, method=grad.method, method.args=list(eps=eps, r=2), ...),silent=TRUE)
                         else  try(gr(pnew, ...),silent=TRUE)
      
      geval <- geval + 1
 
      if (class(gnew)=="try-error" | any(is.nan(gnew)) ){
        lsflag <- 3
        break
        }
 
      s <- pnew - par
      y <- gnew - g
      sts <- sum(s*s)
      yty <- sum(y*y)
      sty <- sum(s*y)
 
      if (method==1) lambda <- min(lmax, max(lmin, sts/sty))
      if (method==2) lambda <- min(lmax, max(lmin, sty/yty))
      if (method==3) lambda <- min(lmax, max(lmin, sqrt(sts/yty)))
 
      if (method==1 & (sts==0 | sty < 0)) lambda <- lmax
      if (method==2 & (sty < 0 | yty == 0)) lambda <- lmax
      if (method==3 & (sts==0 | yty == 0)) lambda <- lmax
 
      par <- pnew
      g   <- gnew
      pg  <- par - g
 
      pg <- if (is.null(project)) try(project.box(pg, lower, upper), silent=TRUE)
  			   else   try(project(pg, ...), silent=TRUE)
 
      if (class(pg) == "try-error" | any(is.nan(pg)) ) {
  	lsflag <- 4
  	break
  	}

      pg <- pg - par
      pg2n <- sqrt(sum(pg*pg))
      pginfn <- max(abs(pg))
 
      f.rep <- (-1)^maximize * f
      if (trace & (iter%%triter == 0))
           cat("iter: ",iter, " f-value: ", f.rep, " pgrad: ",pginfn, "\n")
 
      if (f < fbest) {
  	fbest <- f
  	pbest <- pnew
  	gbest <- pginfn
  	}
 
      }   # while condition loop concludes
 
  if (lsflag==0) {
    if (pginfn <= gtol) conv <- list(type=0, message="Successful convergence")
    if (iter >= maxit)  conv <- list(type=1, message="Maximum number of iterations exceeded")
    } else {
      par <- pbest
      f.rep <- f <- (-1)^maximize * fbest
      pginfn <- gbest
      if (lsflag==1) conv <- list(type=3, message="Failure:  Error in function evaluation")
      if (lsflag==2) conv <- list(type=2, message="Maximum function evals exceeded")
      if (lsflag==3) conv <- list(type=4, message="Failure:  Error in gradient evaluation")
      if (lsflag==4) conv <- list(type=5, message="Failure:  Error in projection")
      }
 
  return(list(par=par, value=f.rep, gradient =pginfn, 
      fn.reduction=(-1)^maximize * (f0 - f), 
      iter=iter, feval=feval, convergence=conv$type, message=conv$message))
  }
 
 
 
 
