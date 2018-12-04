Rcgminu <- function(par, fn, gr, control = list(), ...) {


  ## An R version of the conjugate gradient minimization
  ## using the Dai-Yuan and Hager-Zhang ideas


  #  This version is for unconstrained functions.
  #
  # Input:
  #  par  = a vector containing the starting point
  # fn = objective function (assumed to be sufficeintly
  #   differentiable)
  #  gr = gradient of objective function
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #  control = list of control parameters. Rcgmin uses:
  #           maxit = a limit on the number of iterations (default 500)
  #           (deprecated) maximize = TRUE to maximize the function (default FALSE)
  #           trace = 0 (default) for no output,
  #                  >0 for output (bigger => more output)
  #           eps=1.0e-7 (default) for use in computing numerical gradient approximations.
  #           dowarn=TRUE by default. Set FALSE to suppress warnings.
  # ?? others
  #
  # Output:
  #    A list with components:
  #
  #     par: The best set of parameters found.
  #
  #   value: The value of 'fn' corresponding to 'par'.
  #
  # counts: A two-element integer vector giving the number of
  #   calls to
  # 'fn' and 'gr' respectively. This excludes those calls
  #   needed
  # to compute the Hessian, if requested, and any calls to
  #   'fn'
  # to compute a finite-difference approximation to the
  #   gradient.
  #
  # convergence: An integer code. '0' indicates successful
  #   convergence.
  #          Error codes are
  #          '0' converged
  # '1' indicates that the function evaluation count
  #   'maxfeval'
  #               was reached.
  #          '2' indicates initial point is infeasible
  #
  # message: A character string giving any additional
  #   information returned
  #          by the optimizer, or 'NULL'.
  #
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  #
  #
  #  Author:  John C Nash
  #  Date:  April 2, 2009; revised July 28, 2009
  #         Major rework November 20, 2018 onwards
  #################################################################
  bvec <- par  # copy the parameter vector
  n <- length(bvec)  # number of elements in par vector
  # control defaults -- idea from spg
  ctrl <- cgdefault(n) # use control defaults
  namc <- names(control)
  if (!all(namc %in% names(ctrl))) 
     if(ctrl$trace > 0) warning("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  ctrl[namc] <- control
  npar<-length(par)
  grNULL <- is.null(gr)
  # print(ctrl)
  #############################################
  if (ctrl$maximize) { #!! NOTE
     warning("Rcgmin no longer supports maximize 111121 -- see documentation")
     msg<-"Rcgmin no longer supports maximize 111121"
     ans <- list(par, NA, c(0, 0), 9999, msg)
     return(ans)
  }
  #############################################
  # gr MUST be provided
  if (is.null(gr)) {  # if gr function is not provided STOP (Rvmmin has definition)
     stop("A gradient calculation (analytic or numerical) MUST be provided for Rcgmin") 
  }
  if ( is.character(gr) ) {
  # Convert string to function call, assuming it is a numerical gradient function
     mygr<-function(par=par, userfn=fn, ...){
       do.call(gr, list(par, userfn, ...))
     }
  } else { mygr<-gr }
  ############# end test gr ####################
  ## Set working parameters (See CNM Alg 22)
  if (ctrl$trace > 0) {
    cat("Rcgminu -- J C Nash 2009 - unconstrained version CG min\n")
    cat("an R implementation of Alg 22 with Yuan/Dai modification\n")
  }
  ig <- 0  # count gradient evaluations
  ifn <- 1  # count function evaluations (we always make 1 try below)
  acctol <- ctrl$acctol # acceptable point tolerance ??possibly change later
  reltest <- ctrl$cgoffset # offset for relative test
  stepredn <- ctrl$cgstepredn # factor for shrinking stepsize
  # for ceps
  cyclimit <- min(2.5 * n, 10 + sqrt(n))  # upper bound on when we restart CG cycle
  # This does not appear to be in Y H Dai & Y Yuan, Annals of
  #   Operations Research 103, 33--47, 2001
  # in Alg 22 pascal, we can set this as user. Do we wish to allow that?
  if (ctrl$tol == 0) ctrl$tol <- npar * (npar * .Machine$double.eps)  # for gradient test.  
  ## Note -- integer overflow if n*n*d.eps
  fargs <- list(...)  # function arguments
  if (ctrl$trace > 2) {
    cat("Extra function arguments:")
    print(fargs)
  }
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  ######## check bounds and masks #############
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
    # At this point, we have full bounds in play
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  # Initial function value 
  if (ctrl$trace > 2) {
    cat("Try function at initial point:")
    print(bvec)
  }
  f <- try(fn(bvec, ...), silent = TRUE)  # Compute the function at initial point.
  if (class(f) == "try-error") {
    msg <- "Initial point is infeasible."
    if (ctrl$trace > 0) cat(msg, "\n")
    ans <- list(par, NA, c(ifn, 0), 2, msg)
    names(ans) <- c("par", "value", "counts", "convergence", "message")
    return(ans)
  }
  if (ctrl$trace > 0) cat("Initial function value=", f, "\n")
  fmin <- f # save the value
  if (ctrl$trace > 2) print(bvec)
  # Start the minimization process
  keepgoing <- TRUE
  msg <- "not finished"  # in case we exit somehow
  oldstep <- ctrl$cgoldstep  #?? Why this choice? == .8
  ####################################################################
  fdiff <- NA  # initially no decrease
  cycle <- 0  # cycle loop counter
  kf <- 0
  while (keepgoing) { # main loop -- must remember to break out of it!!
    # (Re)start with steepest descent by zeroing last gradient and last searchdir
    t <- as.vector(rep(0, n))  # zero step vector
    c <- t  # zero 'last' gradient
    while (keepgoing && (cycle < cyclimit)) { ## cycle loop
      cycle <- cycle + 1
      if (ctrl$trace > 0) cat("f, g, kf, cycle:",ifn, " ", ig, " ",kf, " ", cycle, " ", fmin,
                "  last decrease=", fdiff, "\n")
      if (ctrl$trace > 2) {print(bvec); cat("\n")}
      if (ifn > ctrl$maxfeval) {
         msg <- paste("Too many function evaluations (> ", ctrl$maxfeval, ") ", sep = "")
         if (ctrl$trace > 0) cat(msg, "\n")
         ans <- list(par, fmin, c(ifn, ig), 1, msg)  # 1 indicates not converged in function limit
         names(ans) <- c("par", "value", "counts", "convergence", "message")
         return(ans)
      }
      par <- bvec # save best parameters
      ig <- ig + 1 # count gradient evaluations
      if (ig > ctrl$maxit) {
        msg <- paste("Too many gradient evaluations (> ", ctrl$maxit, ") ", sep = "")
        if (ctrl$trace > 0) cat(msg, "\n")
        ans <- list(par, fmin, c(ifn, ig), 1, msg)  # 1 indicates not converged in function or gradient limit
        names(ans) <- c("par", "value", "counts", "convergence", "message")
        return(ans)
      }
      g <- mygr(bvec, ...) # Note: not checked by try(). May want to do that
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
      g1 <- sum(g * (g - c))  # gradient * grad-difference
      g2 <- sum(t * (g - c))  # oldsearch * grad-difference
      gradsqr <- sum(g * g)
      if (ctrl$trace > 1) cat("Gradsqr = ",gradsqr," g1, g2 ",g1," ",g2," fmin=",fmin,"\n")
      c <- g  # save last gradient
      g3 <- 1  # Default to 1 to ensure it is defined -- t==0 on first cycle
      if (gradsqr > ctrl$tol * (abs(fmin) + reltest)) { # ensure we haven't got a "small" gradient
        if (g2 > 0) {
           betaDY <- gradsqr/g2
           betaHS <- g1/g2
           g3 <- max(0, min(betaHS, betaDY))  # g3 is our new 'beta' !! Dai/Yuan 2001, (4.2)
        } #?? What if g2 <= 0?
      } else {
         msg <- paste("Very small gradient -- gradsqr =", gradsqr, sep = " ")
         if (ctrl$trace > 0) cat(msg, "\n")
         keepgoing <- FALSE  # done loops -- break
         break  # to leave inner loop
      }
      if (ctrl$trace > 2) cat("Betak = g3 = ", g3, "\n")
      if (g3 == 0 || cycle >= cyclimit) { # we are resetting to gradient in this case
         if (ctrl$trace > 0) {
           if (cycle < cyclimit) cat("Yuan/Dai cycle reset\n")
             else cat("Cycle limit reached -- reset\n")
         }
         fdiff <- NA
         cycle <- 0
         break  # to quit inner loop
      } else { # drop through if not Yuan/Dai or cycle limit reset
        t <- t * g3 - g  # t starts at zero, later is search vector
        gradproj <- sum(t * g)  # gradient projection
        if (ctrl$trace > 1) cat("Gradproj =", gradproj, "\n")
        accpoint <- FALSE
#        maxstep <- ctrl$steplenn # ?? clean up Do we need?
        if (gradproj < 0) {
  #
  #
  #
  #
  #
  #
  #
          ########################################################
          ####                  Line search                   ####
#          stl <- min(oldstep*ctrl$stinflate, ctrl$cgstep0max)
          if (ctrl$trace > 2) cat("backtracklsq  oldstep=", oldstep,"  ")
          stl <- min(oldstep*ctrl$cgstinflate, ctrl$cgstep0max)
#          stl <- min(oldstep*1.5, 1)
          if (ctrl$trace > 2) cat("stl=", stl,"  ")
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
          kf <- 0
          if (ctrl$trace > 3) cat("qiilev =",ctrl$qiilev,"\n")
          while (! isTRUE(accpoint)) { 
            bvec <- par + stl * t
            if (identical(bvec, par)) {
              msg <- "linesearch: No progress"
              if (ctrl$trace > 0) cat(msg,"\n") 
              stl <- 0 # to ensure flagged
              break
            }
            f <- fn(bvec, ...) # Not testing in try()
            kf <- kf + 1 
            ifn <- ifn + 1 # should test but ...
            if (ctrl$trace > 2) cat("kf=",kf,"   f(",stl,")=",f)
            if (f < fmin + ctrl$qiilev*(abs(fmin)+ctrl$epstol) ) { # try quad point -- note condition
              aa <- (f - fmin - gradproj*stl)/(stl*stl)
              sq <- -gradproj/(2*aa)
              if (ctrl$trace > 2) cat("    aa, sq:",aa,sq,"   ")
              bq <- par + sq*t
              if (! identical(bq, par)) {
                fq <- fn(bq, ...)
                kf <- kf + 1
                ifn <- ifn + 1
                if (ctrl$trace > 2) cat("fq=",fq,"\n")
                if (fq < f) {
                  accpointq <- (fq <= fmin + gradproj * sq * acctol) 
                  if (ctrl$trace > 2) cat("accpointq=",accpointq,"\n")  
                  if (accpointq) {
                    f <- fq
                    stl <- sq
                    bvec <- bq
                  } # end accpointq
                } # end fq < f
              } # end non-identical bq, par
            } # end try quadpoint
            accpoint <- (f <= fmin + gradproj * stl * acctol)
            if (ctrl$trace > 2) cat("   accpoint=", accpoint,"\n")
            if (! isTRUE(accpoint)) stl <- stl*stepredn # backtrack
          } # end backtrack loop
          oldstep <- stl
          fdiff <- fmin - f
          fmin <- f
          par <- bvec
          if (ctrl$trace > 2) cat("Have fmin=",fmin,"  OK accpoint with stl=",stl,"\n")
        } else {
          msg <- "Uphill search direction"
          if (ctrl$trace > 0) cat(msg,"\n") 
          stl <- -1 # flag uphill (probably don't need!!)
        } # gradproj test
        if (stl <= 0) { # not changed on step redn, or uphill
          if (cycle == 1) {
            msg <- " Converged -- no progress on new CG cycle"
            if (ctrl$trace > 0) cat("\n", msg, "\n")
            keekpgoing <- FALSE
            break  #!!
          }
          cycle <- 0 # restart
        }  # end stl <= 0
      }  # end of test on Yuan/Dai condition
#      oldstep <- oldstep*1.5 # ??
#      oldstep <- stl # ?? may want to adjust
      if (ctrl$trace > 2) cat("oldstep=", oldstep,"\n")
      if (oldstep > ctrl$steplenn/1.5) { oldstep <- ctrl$steplenn/1.5}
      if (oldstep < acctol/1.5) { oldstep <- acctol/1.5} #   steplength
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
      if (ctrl$trace > 2) {
        cat("end inner loop -- stl, oldstep = ",stl,oldstep,"\n")
        cat("fmin =",fmin,"\n")
        print(par)
      }
  #   tmp <- readline("again")
    }  # end of inner loop (cycle)        
    if (ctrl$trace > 2) cat("End outer loop, cycle =", cycle, "\n")
  }  # end of outer loop
  msg <- "Rcgminu seems to have converged"
  if (ctrl$trace > 0) cat(msg, "after ",ifn," fn and ",ig," gr evals\n")
  #  par: The best set of parameters found.
  #  value: The value of 'fn' corresponding to 'par'.
  #  counts: number of calls to 'fn' and 'gr' (2 elements)
  # convergence: An integer code. '0' indicates successful
  #   convergence.
  #  message: A character string or 'NULL'.
  #
  #
  #
  ans <- list(par, fmin, c(ifn, ig), 0, msg)
  names(ans) <- c("par", "value", "counts", "convergence", "message")
  return(ans)
}  ## end of Rcgminu
