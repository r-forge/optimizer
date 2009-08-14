Rcgminu <- function( par, fn, gr=NULL, control=list(), ...) {
## An R version of the conjugate gradient minimization
## using the Dai-Yuan ideas
#  This version is for unconstrained functions.
#
# Input:
#  par = a vector containing the starting point
#  fn = objective function (assumed to be sufficeintly differentiable)
#  gr = gradient of objective function
#  control = list of control parameters
#           maxfeval = a limit on the function evaluations (default 10000)
#           maximize = TRUE to maximize the function (default FALSE)
#           trace = 0 (default) for no output, 
#                  >0 for output (bigger => more output)
#
# Output:
#    A list with components: 
#
#     par: The best set of parameters found.
#
#   value: The value of 'fn' corresponding to 'par'.
#
#  counts: A two-element integer vector giving the number of calls to
#          'fn' and 'gr' respectively. This excludes those calls needed
#          to compute the Hessian, if requested, and any calls to 'fn'
#          to compute a finite-difference approximation to the gradient.
#
# convergence: An integer code. '0' indicates successful convergence.
#          Error codes are
#          '0' converged 
#          '1' indicates that the function evaluation count 'maxfeval'
#               was reached.
#          '2' indicates initial point is infeasible
#
# message: A character string giving any additional information returned
#          by the optimizer, or 'NULL'.
#
#  Author:  John C Nash
#  Date:  April 2, 2009; revised July 28, 2009
#################################################################
  # control defaults -- taken from spg -- need to adjust for Rcgmin and Rvmmin
  ctrl<-list(maxfeval=10000, maximize=FALSE, trace=0)
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])     
  ctrl[namc ] <- control
  maxfeval <- ctrl$maxfeval  # limit on function evaluations
  maximize <- ctrl$maximize  # TRUE to maximize the function
  trace<- ctrl$trace         # 0 for no output, >0 for output (bigger => more output)
  fargs <- list(...)         # the ... arguments that are extra function / gradient data
#############################################
# local function defined only when user does not specify a gr
# Simple gr numerical approximation. Using fn, fmin and eps from calling env.  	
  if (is.null(gr)) {
       warning("WARNING: Numerical gradients are inappropriate for Rcgminu")
       eps<-1.0e-7 # should we use a different value?
       gr <- function(par, ...) {
    	df <- rep(NA,length(par))
        teps<-eps*(abs(par)+eps)
    	for (i in 1:length(par)) {
    	  dx <- par
    	  dx[i] <- dx[i] + teps[i]
    	  df[i] <- (fn(dx, ...) - fmin)/teps[i]
    	 }
    	df
	}
   }
############# end gr ########################
## Set working parameters (See CNM Alg 22)
  if (trace>0) cat("Rcgminu -- J C Nash 2009 - an R implementation of Alg 22 in Yuan/Dai version\n")
  bvec <- par # copy the parameter vector
  n <- length(bvec) # number of elements in par vector
  ig <- 0 # count gradient evaluations
  ifn <- 1 # count function evaluations (we always make 1 try below)
  stepredn <- 0.15 # Step reduction in line search
  acctol <- 0.0001 # acceptable point tolerance
  reltest <- 100.0 # relative equality test
  setstep <- 1.75 # step increase factor
  accpoint <- as.logical(FALSE) # so far do not have an acceptable point
  fail <- as.logical(FALSE) # Method hasn't yet failed on us!
  cyclimit <- min(n,10+sqrt(n)) # upper bound on when we restart CG cycle 
  # This does not appear to be in Y H Dai & Y Yuan, Annals of Operations Research 103, 33–47, 2001 aor01.pdf
  intol <- .Machine$double.eps # machine precision from .Machine variable 
  ## in Alg 22 pascal, we can set this as user
  tol <- n*(n*intol) # try different Note -- integer overflow if n*n*intol
#  tol <- n * intol * sqrt(intol) # a tolerance for gradient test
  f <- try(do.call("fn", append(list(bvec), fargs )), silent=TRUE) # Compute the function.
  if ( class(f) == "try-error") { 
     msg<-"Initial point is infeasible."
     if(trace > 0) cat(msg,"\n")
     ans<-list(x, NA, c(ifn, 0), 2, msg)
     return(ans)
  } 
  fmin<-f
  keepgoing<-TRUE
  notconv <-TRUE
  msg<-"not finished" # in case we exit somehow
  oldstep<-1.0
  while (keepgoing && notconv) {  ## main loop -- must remember to break out of it!
    t<-as.vector(rep(0,n)) # zero step vector
    c<-t # zero "last" gradient
    cycle <- 0 # cycle loop counter
    count<-0
    while (keepgoing && (cycle < cyclimit)) { ## cycle loop
      cycle<-cycle+1
      if (trace > 0) cat(ifn," ",ig," ",cycle," ",fmin)
      if( trace>2) {
        print(bvec)
        cat("\n")
      }
      if (ifn > maxfeval) {
        msg<-paste("Too many function evaluations (> ",maxfeval,") ", sep='')
        if (trace > 0 ) cat(msg,"\n")
        ans<-list(x, fmin, c(ifn, ig), 1, msg) # 1 indicates not converged in function limit
        names(ans)<-c("par", "value", "counts", "convergence", "message")
        return(ans)
      }
      x<-bvec # save best parameters
      ig<-ig+1
      g<-gr(bvec, ...)
      g1<-sum(g*(g-c)) # gradient * grad-difference
      g2<-sum(t*(g-c)) # oldsearch * grad-difference
      gradsqr<-sum(g*g)
      if (trace > 1) cat("Gradsqr = ", gradsqr," g1, g2 ", g1, " ", g2," fmin=",fmin,"\n")
      c<-g # save last gradient
      g3<-1.0 # Default to 1 to ensure it is defined
      if (gradsqr > tol*(abs(fmin)+reltest) ) {
         if (g2 > 0.0) {
            betaDY<-gradsqr/g2
            betaHS<-g1/g2
            g3 <-max( 0.0, min(betaHS, betaDY) ) # g3 is our new 'beta'
         } 
      } else {
	 msg<-paste("Very small gradient -- gradsqr =",gradsqr,sep=' ')
         if (trace > 0) cat(msg,"\n")
         keepgoing<-FALSE
         break # to leave inner loop
      }    
      if(trace > 2) cat("Betak = g3 = ",g3,"\n")
      if (g3 == 0.0) {  # we are resetting to gradient in this case 
          if (trace > 0) cat("Yuan/Dai cycle reset\n")
          cycle<-cyclimit+1
          break # to quit inner loop
      } else {
        # drop through if not Yuan/Dai cycle reset
        t <- t * g3 - g # t starts at zero, later is step vector
        gradproj <- sum(t*g) # gradient projection
        if(trace > 1) cat("Gradproj =",gradproj,"\n")
########################################################
####                  Line search                   ####
        OKpoint<-FALSE
        steplength <- oldstep
        f <- fmin
        changed<-TRUE
        while ((f >= fmin) && changed) {
          bvec<-x+steplength*t
          changed <- ( !  identical((bvec+reltest), (x+reltest)) )
          if( changed ) { # compute newstep, if possible
            f <- fn(bvec, ...) # Because we need the value for linesearch, don't use try()
            #  instead preferring to fail out, which will hopefully be unlikely.
            ifn<-ifn+1
            if (f < fmin) {
              f1 <- f # ?? save just in case
            } else {
              steplength<-steplength*stepredn
              if (trace > 0) cat("*")
            }
          }
        }  # end while
        changed1<-changed
        if (changed1) {
          ## ?? should we check for reduction? or is this done in if (newstep >0) ?
          newstep <- 2.0*(f - fmin - gradproj * steplength) # JN 081219 change
          if (newstep > 0) {
              newstep = -(gradproj * steplength * steplength / newstep) 
          }
          bvec<-x+newstep*t
          changed <- ( ! identical((bvec+reltest), (x+reltest)) )
          if (changed) {
              f <- fn(bvec, ...)
              ifn<-ifn+1
          }
          if(trace > 1) cat("fmin, f1, f: ",fmin, f1, f,"\n")
          if (f < min(fmin, f1)) {
             # success
             OKpoint<-TRUE
             accpoint <- ( f <= fmin + gradproj * newstep * acctol) 
             fmin<-f
          } else { 
            if (f1 < fmin) {
	        bvec<-x+steplength*t # reset best point
                accpoint <- ( f1 <= fmin + gradproj * steplength * acctol) 
                OKpoint<-TRUE ## may need to fix up
                fmin<-f1
            } else {
                # no reduction
                accpoint<-FALSE
            } # f1<?fmin
          } # f < min(f1, fmin)
          if(trace > 1) cat("accpoint = ",accpoint," OKpoint = ",OKpoint,"\n")
          if (! accpoint) {
              msg<-"No acceptable point -- exit loop"
              if (trace > 0) cat("\n",msg,"\n")
              keepgoing<-FALSE
             # break
          }
        } # changed1
        else { # not changed on step redn
          if (cycle == 1) { 
            msg<-" Converged -- no progress on new CG cycle"
            if (trace > 0) cat("\n",msg,"\n")
            keekpgoing<-FALSE
            notconv<-FALSE # mark as converged
          } 
        } # end else
      } # end of test on Yuan/Dai condition      
    } # end of inner loop (cycle)
    if (trace > 1) cat("Break inner loop, cycle =",cycle,"\n")
#     temp<-readline("continue?")
    if (oldstep < acctol) { oldstep <- acctol } 
    oldstep = setstep * newstep # change to newstep fro steplength
    if (oldstep > 1.0) { oldstep <- 1.0 }
  } # end of outer loop
  if (trace>0) cat("Seem to be finished Rcgminu\n\n\n")
 #  par: The best set of parameters found.
 #  value: The value of 'fn' corresponding to 'par'.
 #  counts: number of calls to 'fn' and 'gr' (2 elements)
 #  convergence: An integer code. '0' indicates successful convergence.
 #  message: A character string or 'NULL'.
  ans<-list(x, fmin, c(ifn, ig), 0, msg)
  names(ans)<-c("par", "value", "counts", "convergence", "message")
  return(ans)
} ## end of Rcgmin

