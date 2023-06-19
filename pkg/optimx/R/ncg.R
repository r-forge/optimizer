ncg <- function(par, fn, gr, lower=NULL, upper=NULL, bdmsk = NULL, control = list(), ...) {
  ## ONLY minimizes -- use efn, egr with scaling already set. Call from optimr
  ## Get approximation into egr in optimr also.
    if (is.null(gr)) stop("A gradient calculation (analytic or numerical) MUST be provided for ncg") 

   tryxtra <- FALSE # CONTROL ?? for trying an extra step
## Feb 15 -- working but inefficient for exrosen
## Revised 2206 to use simpler linesearch
    ## An R version of the conjugate gradient minimization using the Dai-Yuan ideas
    # ncg.R 20220212 JN
    npar<-length(par)
    ctrl <- ctrldefault(npar)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$tol == 0) tol <- npar * (npar * .Machine$double.eps)  
       # for gradient test.  Note -- integer overflow if n*n*d.eps
    else tol<-ctrl$tol
    maxit <- ctrl$maxit  # limit on function evaluations
    trace <- ctrl$trace  # 0 for no output, >0 for output (bigger => more output)
    stepredn <- ctrl$stepredn
    eps <- ctrl$eps
    dowarn <- ctrl$dowarn  #
    bounds <- ! all(as.logical(bdmsk))
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) cat("ncg -- J C Nash 2022 - bounds constrained version of CG\n")
    bvec <- par  # copy the parameter vector
    n <- length(bvec)  # number of elements in par vector
    maxfeval <- round(sqrt(n + 1) * maxit)  # change 091219
    ig <- 0  # count gradient evaluations
    ifn <- 1  # count function evaluations (we always make 1 try below)
    acctol <- ctrl$acctol  # acceptable point tolerance
    reltest <- ctrl$reltest  # relative equality test 
    ceps <- .Machine$double.eps * reltest
    pceps <- ceps*max(abs(bvec))
    accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
    cyclimit <- min(2.5 * n, 10 + sqrt(n))  #!! upper bound on when we restart CG cycle
    # set default masks if not defined. Normally set by optimr.
    if (is.null(bdmsk)) { bdmsk <- rep(1, n) }
    if (trace > 2) { cat("bdmsk:"); print(bdmsk)}
    # Initial function value -- may NOT be at initial point specified by user.
    if (trace > 2) {cat("Try function at initial point:"); print(bvec) }
    f <- try(fn(bvec, ...), silent = TRUE)  # Compute the function at initial point.
    if (trace > 0) {cat("Initial function value=", f, "\n") }
    if (inherits(f,"try-error")) {
      msg <- "Initial point is infeasible."
      if (trace > 0) cat(msg, "\n")
      ans <- list(par, NA, c(ifn, 0), 2, msg, bdmsk)
      names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
      return(ans)
    }
    fmin <- f
    if (trace > 0) cat("Initial fn=", f, "\n")
    if (trace > 1) print(bvec)
    # Start the minimization process
    keepgoing <- TRUE
    msg <- "not finished"  # in case we exit somehow
    oldstep <- 0.8  #!! 2/3 #!!? WHY? formerly 0.8
    ####################################################################
    fdiff <- NA  # initially no decrease
    cycle <- 0  # !! cycle loop counter
    haveg <- FALSE
    while (keepgoing) {
      # main loop -- must remember to break out of it!!
      t <- as.vector(rep(0, n))  # zero step vector
      c <- t  # zero 'last' gradient
      while (keepgoing && (cycle < cyclimit)) {
        ## cycle loop
        cycle <- cycle + 1
        if (trace > 0) cat(ifn, " ", ig, " ", cycle, " ", fmin, "  last decrease=", fdiff, "\n")
        if (trace > 1) { print(bvec) }
        if (! haveg) {
          par <- bvec  # save best parameters
          ig <- ig + 1
          if (ig > maxit) {
            msg <- paste("Too many gradient evaluations (> ", maxit, ") ", sep = "")
            if (trace > 0) cat(msg, "\n")
            ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  
            # 1 indicates not converged in function or gradient limit
            names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
            return(ans)
          }
          g <- gr(bvec, ...) # gradient
          haveg <- TRUE
        }
        if (bounds) { ## Bounds and masks adjustment of gradient
          ## first try with looping -- later try to vectorize
          if (trace > 2) { cat("bdmsk:"); print(bdmsk) }
          for (i in 1:n) {
            if ((bdmsk[i] == 0)) { g[i] <- 0 } # masked, so gradient component is zero
            else {
              if (bdmsk[i] != 1) {
                if ((bdmsk[i] + 2) * g[i] < 0) { # test for -ve gradient at upper bound, +ve at lower bound
                  g[i] <- 0  # active mask or constraint and zero gradient component
                }
                else {
                  bdmsk[i] <- 1  # freeing parameter i
                  if (trace > 1) cat("freeing parameter ", i, "\n")
                }
              }
            }
          }  # end masking loop on i
          if (trace > 2) {cat("bdmsk adj:\n"); print(bdmsk); cat("proj-g:\n"); print(g) }
        }  # end if bounds
        ## end bounds and masks adjustment of gradient
        g1 <- sum(g * (g - c))  # gradient * grad-difference
        g2 <- sum(t * (g - c))  # oldsearch * grad-difference
        gradsqr <- sum(g * g)
        if (trace > 1) { cat("Gradsqr = ", gradsqr, " g1, g2 ", g1, " ", g2, " fmin=", fmin, "\n") }
        c <- g  # save last gradient
        g3 <- 1  # !! Default to 1 to ensure it is defined -- t==0 on first cycle
        if (gradsqr > tol * (abs(fmin) + reltest) ) {
          if (g2 > 0) {
            betaDY <- gradsqr/g2
            betaHS <- g1/g2
            g3 <- max(0, min(betaHS, betaDY))  # g3 is our new 'beta' !! Dai/Yuan 2001, (4.2)
          }
        }
        else {
          msg <- paste("Very small gradient -- gradsqr =", gradsqr, sep = " ")
          if (trace > 0) cat(msg, "\n")
          keepgoing <- FALSE  # done loops -- should we break?
          break  # to leave inner loop
        }
        if (trace > 2) cat("Betak = g3 = ", g3, "\n")
        if (g3 == 0 || cycle >= cyclimit) { # we are resetting to gradient in this case
          if (trace > 0) {
            if (cycle < cyclimit) cat("Yuan/Dai cycle early reset\n")
            else cat("Cycle limit reached or forced -- reset\n")
          }
          fdiff <- NA
          cycle <- 0 # but haveg == TRUE
          break  #!!
        }
        else { # drop through if not Yuan/Dai cycle reset
          gpoldt <- sum(g * t)# ?? testcode
          if (cycle > 1) {
              if (trace > 1) cat(cycle, "last gradproj=",gradproj,"  new g old t=",gpoldt,"\n")# ?? testcode
          }
          t <- t * g3 - g  # t starts at zero, later is step vector
          gradproj <- sum(t * g)  # gradient projection
          if (trace > 1) cat("Gradproj =", gradproj, "\n")
          if (bounds) { ## Adjust search direction for masks
            if (trace > 2) { cat("t:\n"); print(t) }
            t[which(bdmsk <= 0)] <- 0  # apply mask constraint
            if (trace > 2) { cat("adj-t:\n"); print(t) }
            ## end adjust search direction for masks
          }  # end if bounds
          # Why do we not check gradproj size??
          accpoint <- FALSE
          f <- fmin # and as large as fmin
          if (trace > 1) cat("Start linesearch with oldstep=", oldstep, "\n")
          stlen <- oldstep * 1.5  #!! try a bit bigger to match Rcgmin
          stepstrt<-stlen
          stepmax <- ctrl$bigval # large number
          if (bounds) { # Box constraint -- adjust step length
            # ?? put in function (possibly at top of this routine to use local vars
            for (i in 1:n) { # loop on parameters -- vectorize?
              if ((bdmsk[i] == 1) && (abs(t[i]) > pceps)) { # only free params and search != 0
                if (t[i] < 0) { # going downhill. Look at lower bound
                  trystep <- (lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
                }
                else { # going uphill, check upper bound
                  trystep <- (upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
                }
                # if (trace > 3) cat("stlen, trystep:", stlen, trystep, "\n")
                stepmax <- min(stepmax, trystep)  # reduce as necessary
              }  # end stlen reduction
              else {t[i] <- 0} # small dirn gets set to zero just in case 
            }  # end loop on i to reduce step length
             # end box constraint max step length
          }  # end if bounds
          bdlim <- (stepmax < stepstrt) # TRUE if bound limits step ?? remove?
          if (trace > 1) cat("reset stepmax (",bdlim,") = ", stepmax, "\n")
          stlen <- min(stepstrt, stepmax) # use limit
          # ?? may not need all these indicators
          changed <- TRUE  # Need to set so loop will start
          accpoint <- FALSE
          lopt <- FALSE # point is NOT acceptable, nor lower
          lstry <- TRUE # line search not finished
          f0 <- fmin # save starting fval
          ## LINE SEARCH
          # 220628 retry -- note stlen already set by 0.8 == oldstep and multipy by 1.5 
          # ?? would 1.25 be better --> 1 as start
          while (lstry) { # keep going until we CANNOT change the point
            bvec <- par + stlen *t
            changed <- any(par+reltest != bvec+reltest)
            if (changed){
              ifn <- ifn + 1
              if (ifn > maxfeval) {
                msg <- paste("Too many function evaluations (> ", maxfeval, ") ", sep = "")
                if (trace > 0) cat(msg, "\n")
                ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  # 1 indicates not converged in function limit
                names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
               return(ans)
              }
              f <- fn(bvec, ...)  # Because we need the value for linesearch, don't use try()??
              if (is.na(f) || (!is.finite(f))) {
                       warning("ncg - undefined function")
                       f <- .Machine$double.xmax
              }
              if (trace > 2) cat(ifn," f=",f," at step ",stlen,"\n")
              accpoint <- (f < f0 + gradproj * stlen * acctol) # < not <=
              if (accpoint) { # have an OK point, but should we try further?
                fmin <- f
                if (tryxtra) {
                  if (f < (f0 - gradproj*sbest) ) { # f < expected from linear extrapolation
                    newstep <- min(stepmax, 2.0*stlen) # try doubling
                  } 
                  else {
                    aa <- (f - f0 - gradproj*sbest)/sbest^2
                    newstep <- -gradproj/(2*aa)
                  }
                  if ((newstep+reltest) != (stlen+reltest)) {
                    changed <- TRUE # just to be save
                    if (trace>1) cat("Extra step size=",newstep)
                    bvec <- par + newstep*t
                    newf <- fn(bvec, ...)
                    ifn <- ifn + 1
                    if (trace > 1) cat(" newfval =",newf,"\n")
                    if (newf < f) { # ?? is it acceptable? Assume yes
                       accnew <- (newf <= f0 + gradproj * newstep * acctol) 
                       if (trace > 1) cat("newf=",newf," acnew =",accnew,"\n")
                       sbest <- newstep
                       fmin <- newf
                    }
                    else {
                        if (trace >1) cat("newstep fails to change step\n")
                    }
                  } # newstep changes prams
                } # tryextra
                lstry <- FALSE # done. bvec should be best points
                if (trace > 2) {cat("new pt:"); print(bvec)}
              } # accpoint
              else { # pst not acceptable, decrease step
                stlen <- stlen*stepredn
                if (trace > 0) cat("*") 
              } # NOT lopt
            }
            else { # no change - par has best parameters, f0 best fn
              if (trace > 1) cat("No change in params in line search\n")
              lstry <- FALSE
#              bvec <- par # ?? is this needed
              # break # out of this while loop to end search
            } # end changed
#            tmp <- readline("end lstry loop")
          } # while lstry
        }
#         #### End line search ####
        if (bounds && changed) { ## Reactivate constraints? -- should check for infinite bounds!!?
          for (i in 1:n) {
            if (bdmsk[i] == 1) { # only interested in free parameters
              if (is.finite(lower[i])) {# JN091020 -- need to use abs in case bounds negative
                if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) { # are we near or lower than lower bd
                  if (trace > 2) cat("(re)activate lower bd ", i, " at ", lower[i], "\n")
                  bdmsk[i] <- -3
                }  # end lower bd reactivate
              }
              if (is.finite(upper[i])) { # JN091020 -- need to use abs in case bounds negative
                if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) { # are we near or above upper bd
                  if (trace > 2) cat("(re)activate upper bd ", i, " at ", upper[i], "\n")
                  bdmsk[i] <- -1
                }  # end upper bd reactivate
              }
            }  # end test on free params
          }  # end reactivate constraints loop
        }  # end if bounds
        oldstep <- stlen
        haveg <- FALSE # to force rerun
        if (oldstep < acctol) { oldstep <- acctol }  
        if (oldstep > 1) { oldstep <- 1 }
        if (lopt && (! accpoint)) { # reduction -- use steepest descent
          if (trace > 1) cat("cycle ",cycle," lower but not acceptable point\n")
          cycle <- 0  # force reset
          haveg <- FALSE # forece recalc
          break # out of inner loop
        }
      }  # end of inner loop (cycle)
      if (trace > 1) cat("End inner loop, cycle =", cycle, "\n")
    }  # end of outer loop
    msg <- "ncg seems to have converged"
    if (trace > 0) 
        cat(msg, "\n")
    #  par: The best set of parameters found.
    #  value: The value of 'fn' corresponding to 'par'.
    #  counts: number of calls to 'fn' and 'gr' (2 elements)
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #  message: A character string or 'NULL'.
    ans <- list(par, fmin, c(ifn, ig), 0, msg, bdmsk)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message", "bdmsk")
    return(ans)
}  ## end of ncg
