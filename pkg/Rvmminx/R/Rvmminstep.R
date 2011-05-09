Rvmminstep <- function( par, fn, gr=NULL, lower=NULL, upper=NULL, bdmsk=NULL, control=list(), ...) {  #1
## ?? return bdmsk with solution? or some other report on bounds??
## An R version of the Nash version of Fletcher's Variable Metric minimization
#    including bounds constraints on parameters.
# This version uses a simple backtracking line search, which if successful, is then 
#    optionally followed by a quadratic minimum if the successful backtrack 
#    point is above the gradient extrapolation line.
# Input:
#  par = a vector containing the starting point
#  fn = objective function (assumed to be sufficeintly differentiable)
#  gr = gradient of objective function
#  lower = vector of lower bounds on parameters
#  upper = vector of upper bounds on parameters
#     Note: free parameters outside bounds will be adjusted to bounds.
#  bdmsk = control vector for bounds and masks. Parameters for which bdmsk are 1
#         are unconstrained or "free", those with bdmsk 0 are masked i.e., fixed.
#         For historical reasons, we use the same array as an indicator that a
#         parameter is at a lower bound (-3) or upper bound (-1)
#     Note: masked parameters out of bounds are NOT adjusted to bounds.
#  control = list of control parameters
#           maxit = a limit on the gradient evaluations (default 1500)
#                  ?? do we want to make it vary with n??
#           maxfeval = a limit on the function evaluations (default 10000)
#           maximize = TRUE to maximize the function (default FALSE)
#           trace = 0 (default) for no output, 
#                  >0 for output (bigger => more output)
#           eps=1.0e-7 (default) for use in computing numerical gradient
#                  approximations.
#           usenumDeriv=FALSE default. TRUE to use numDeriv for numerical
#                  gradient approximations.) 
#           gradtol = tolerance for small gradient projection. The algorithm
#                  will terminate if the gradient projection magnitude is 
#                  less than or equal to gradtol
#           acctol = tolerance to use for acceptable point test
#           reltest = offset used for relative test of equality
#           stepredn = value to use for reducing the stepsize in backtrack 
#                  search for lower point
#           tryqmin = TRUE (default) if we want to try the quadratic minimization
#                  after backtrack line search, FALSE to skip this attempt.
#           dowarn=TRUE by default. Set FALSE to suppress warnings.
#
##
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
# termination: An integer code. '0' indicates successful termination.
#          Other codes are
#          '1' indicates that the iteration limit 'maxit' or the function
#               evaluation limit maxfeval have been reached.
#          ?? do we want to indicate small gradproj
#          '20' indicates inadmissible input parameters for function (bounds may be OK)
#          '21' indicates inadmissible intermediate parameters for function
#
# message: A character string giving any additional information returned
#          by the optimizer, or 'NULL'.
#
# bdmsk:   Returned index describing the status of bounds and masks at the
#          proposed solution. Parameters for which bdmsk are 1 are unconstrained
#          or "free", those with bdmsk 0 are masked i.e., fixed. For historical
#          reasons, we indicate a parameter is at a lower bound using -3 
#          or upper bound using -1.
#
#  Author:  John C Nash
#  Date:  April 3, 2009
#
# Things done:
# 20110424 added qmin to search
#################################################################
  # control defaults -- some taken from spg -- ?? need to adjust for Rcgmin and Rvmmin
  ctrl <- list( maxit=50, maxfeval=300, maximize=FALSE, trace=0, 
           eps=1.0e-7, usenumDeriv=FALSE, acctol=0.0001, 
           reltest=1000.0, stepredn=0.2, gradtol=-1, tryqmin=TRUE, dowarn=TRUE)
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  if ((!is.null(control$trace)) && (control$trace>0)){ #
     cat("Names in control:")
     print(namc)
  }
  tmp<-readline("continue")     
  ctrl[namc] <- control      # 
  maxit    <- ctrl$maxit     #
  maxfeval <- ctrl$maxfeval  #
  maximize <- ctrl$maximize  # TRUE to maximize the function
  trace    <- ctrl$trace     #
  eps      <- ctrl$eps       #
  reltest  <- ctrl$reltest # relative equality test  
  stepredn <- ctrl$stepredn # Step reduction in line search
  acctol <- ctrl$acctol # acceptable point tolerance
  gradtol <- ctrl$gradtol # gradient projection tolerance
  tryqmin <- ctrl$tryqmin # try the quadratic min
  usenumDeriv<-ctrl$usenumDeriv # to force use of numDeriv package
  dowarn   <- ctrl$dowarn    #
  grNULL <- is.null(gr)      # if gr function is not provided, we want to use numDeriv
  fargs <- list(...)         # the ... arguments that are extra function / gradient data
  msg <- "" # ensure msg is set, even to a null string

#################################################################
## Set working parameters (See CNM Alg 22)
  if (trace > 0) cat("Rvmmin -- J C Nash 2009 - an R implementation of Alg 21\n")
  bvec <- par # copy the parameter vector
  n <- length(bvec) # number of elements in par vector
  if(trace>0) {
	cat("Problem of size n=",n,"  Extra arguments:\n")
	print(fargs)
  }
  ifn <- 1 # count function evaluations
  ceps<-.Machine$double.eps*reltest # only used for test of nearness to bounds for re-est.
#############################################
# local function defined in order to deal with out of bounds functions/parameters
# ?? add exceeding function count inside and change attributes??
   myfn <- function(par, ...) {
      # cat("In myfn, maximize=",maximize," par =")
      # print(par)
      # cat("list(par):")
      # print(list(par))
      testf<-try(tryf <- fn(par, ...), silent=TRUE) 
      # cat("testf=")
      # print(testf)
      # Compute the function.
      if ( (class(testf) == "try-error") | is.na(tryf) | is.null(tryf) | is.infinite(tryf)) { 
         tryf<- if (maximize) { -.Machine$double.xmax } else { .Machine$double.xmax }
         attr(tryf,"inadmissible")<-TRUE
      } else {
         attr(tryf,"inadmissible")<-FALSE
      }
      # cat("In myfn, function attribute is ",attr(tryf,"inadmissible"),"\n")
      if (maximize) -tryf else tryf # handle the maximization
   }
############# end myfn ##########################
  if (is.null(gr)) {# local function defined only when user does not specify a gr
       if (dowarn) warning("Numerical gradients may be inappropriate for Rvmmin")
       if ( usenumDeriv ) { # external gradient
         require("numDeriv")
         mygr<-function(par, maximize=FALSE, ...) {
           gg<-grad(myfn, par, ...)
           gg
         } # end local gr
       } else { # using local numerical gradient
         # Simple gr numerical approximation. Using fn, fmin and eps from calling env.  	
         # cat("Set gr to simple function with eps=",eps,"\n")
         mygr <- function(par, maximize=FALSE, ...) {
           fbase<-myfn(par,...) # ensure we have right value, may not be necessary
           df <- rep(NA,length(par))
           teps<-eps*(abs(par)+eps)
           for (i in 1:length(par)) {
             dx <- par
             dx[i] <- dx[i] + teps[i]
             df[i] <- (myfn(dx, ...) - fbase)/teps[i]
           }
           df
         } # end mygr simple
       } # end else simple numgrad
  } else {  # not is.null(gr)
       mygr<-function(par, maximize=FALSE, ...) {
           df<-gr(par,...)
           if (maximize) -df else df
       }
  } # end redefine gradient
############## end gradient redefine ###############
# set default masks if not defined
  if (is.null(bdmsk)) {
       bdmsk<-rep(1,n)
  }
  if (trace > 3) {
     cat("bdmsk:")
     print(bdmsk)
  }
# check if there are bounds
  if(is.null(lower) || ! any(is.finite(lower))) nolower=TRUE else nolower=FALSE
  if(is.null(upper) || ! any(is.finite(upper))) noupper=TRUE else noupper=FALSE
  if(nolower && noupper) bounds=FALSE else bounds=TRUE
  if(any(bdmsk == 0)) masks=TRUE else masks=FALSE
  if (trace > 2) {
     cat("Bounds: nolower = ",nolower,"  noupper = ",noupper," bounds = ",bounds,"\n")
     cat("Masks: ",masks,"\n")
  }
  if(nolower) lower<-rep(-Inf,n)
  if(noupper) upper<-rep(Inf,n)
######## check bounds and masks #############
## NOTE: do this inline to avoid call (??should we change this?)
  if (bounds) {
  ## tmp<-readline("There are bounds ")
# Make sure to expand lower and upper
   if(! nolower & (length(lower)<n)) {
        ## tmp<-readline("Check length lower ")
        if (length(lower)==1) { lower<-rep(lower,n) } else { stop("1<length(lower)<n") }
   } # else lower OK
   if(! noupper & (length(upper)<n)) {
        ## tmp<-readline("Check length upper ")
        if (length(upper)==1) { upper<-rep(upper,n) } else { stop("1<length(upper)<n") }
   } # else upper OK
# At this point, we have full bounds in play
# This implementation as a loop, but try later to vectorize
   for (i in 1:n) {
#       cat("i = ",i,"\n")
       if (bdmsk[i] == 0) { # NOTE: we do not change masked parameters, even if out of bounds
           ## tmp<-readline("Masked parameter ")
           if(! nolower) {
              if(bvec[i]<lower[i]) { 
                wmsg<-paste(bvec[i]," = MASKED x[",i,"] < lower bound = ",lower[i],sep='') 
                if (dowarn) warning(wmsg)
              }
           }
           if(! noupper) {
              if (bvec[i]>upper[i]) { 
                wmsg<-paste(bvec[i]," = MASKED x[",i,"] > upper bound = ",upper[i],sep='') 
                if (dowarn) warning(wmsg)
              }
           }
       } else { # not masked, so must be free or active constraint
           ## tmp<-readline(" Not masked parameter ")
           if(! nolower){
              if (bvec[i]<=lower[i]) { # changed 090814 to ensure bdmsk is set
                wmsg<-paste("x[",i,"], set ",bvec[i]," to lower bound = ",lower[i],sep='')
                if (dowarn) warning(wmsg)
                bvec[i]<-lower[i]
                bdmsk[i] <- -3 # active lower bound
              }
           }
           if(! noupper) {
              if(bvec[i]>=upper[i]) { # changed 090814 to ensure bdmsk is set 
                wmsg<-paste("x[",i,"], set ",bvec[i]," to upper bound = ",upper[i],sep='')
                if (dowarn) warning(wmsg)
                bvec[i]<-upper[i]
                bdmsk[i] <- -1 # active upper bound
              }
           }
       } # end not masked
    } # end loop for bound/mask check
  }
############## end bounds check #############
  f <- myfn(bvec, ...)
  if (attr(f,"inadmissible")) {
     msg<-"Initial point gives inadmissible function value"
     termination<-20
     if(trace>0) cat(msg,"\n")
     ans<-list(bvec, NA, c(ifn, 0), 0, termination, msg, bdmsk) # 
     names(ans)<-c("par", "value", "counts", "termination", "message","bdmsk")
     return(ans)
  }
  fstart<-f
  if (gradtol < 0) { gradtol <- sqrt(.Machine$double.eps)*abs(fstart) }
  if (trace>1) cat("gradtol to use is ",gradtol,"\n")
  if (trace>0) cat("Initial fn=",f,"\n")
  if (trace>2) print(bvec)
  keepgoing<-TRUE # to ensure loop continues until we are finished
  ig <- 1 # count gradient evaluations
  ilast<-ig # last time we used gradient as search direction
  fmin<-f # needed for numerical gradients
  cat("check on gradient at bvec:")
  print(bvec)
  print(mygr)
  cat("function is ",f,"\n")
  g<-mygr(bvec, ...) # Do we need to use try() ?? Possible not
  if (trace>2) {
    cat("g:")
    print(g)
  }
  oldstep<-1.0
  termination<- -1 # to indicate NOT having good termination
  while (keepgoing) {  ## main loop -- must remember to break out of it!
    if (ilast == ig) { # reset the approx. inverse hessian B to unit matrix 
      B<-diag(1,n,n) # create unit matrix of order n
      if (trace>2) cat("Reset Inv. Hessian approx at ilast = ",ilast,"\n")
    } 
    fmin<-f
    if (trace>0) cat(" ",ifn," ",ig," ",fmin,"\n")
    par<-bvec # save parameters
    c<-g # save gradient
    ## Bounds and masks adjustment of gradient ##
    ## first try with looping -- later try to vectorize
    if (bounds || masks) {
      if (trace>2) { 
         cat("bdmsk:")
         print(bdmsk)
      }
      for (i in 1:n) {
         if( (bdmsk[i]==0) ) {
            g[i]<-0
         } else {
            if (bdmsk[i]==1) {
               if (trace>2) cat("Parameter ",i," is free\n")
            } else {
               if ( (bdmsk[i]+2)*g[i]<0.0 ) {
                  g[i] <- 0. # active mask or constraint
               } else {
                  bdmsk[i] <- 1 # freeing parameter i
                  if (trace>1) cat("freeing parameter ",i,"\n")
               }
            }
         }
      } # end masking loop on i
      if (trace > 2) {
         cat("bdmsk adj:")
         print(bdmsk)
         cat("proj-g:")
         print(g)
      }
      ## end bounds and masks adjustment of gradient
    } # if bounds
    t<-as.vector(-B%*%g) # compute search direction
    if (trace > 2) {
      cat("t:")
      print(t)
    }
    t[which(bdmsk<=0)]<-0 # apply masks and box constraints
    if (trace > 2) {
       cat("adj-t:")
       print(t)
    }
    gradproj <- sum(t*g) # gradient projection
    if (trace>1) cat("Gradproj =",gradproj,"\n")
    accpoint<-FALSE # Need this BEFORE gradproj test
    if (gradproj < 0) { # Must be going downhill
#     compute maxstep -- inlined from bmchk package
      maxstep<- +Inf
      if (bounds) {
         sn<-rep(+Inf, n)
         sp<-sn
         # these are the max steps to bounds
	 # Note that we need to recognize idxn represents negative t
         idxn<-which(t<0)
         idxp<-which(t>0)
         cat("idxn, idxp:")
         print(idxn)
         print(idxp)
         suppressWarnings(sn[idxn]<-((lower-par)/t)[idxn])
         suppressWarnings(sp[idxp]<-((upper-par)/t)[idxp]) # Awkward!
	 cat("sn, sp:")
	 snsp<-c(sn,sp)
	 snsp<-snsp[which(snsp>0)]
	 print(sn)
         print(sp)
         maxstep<-min(snsp)
      }
      if (trace>2) cat("maxstep =",maxstep,"\n")
      if (abs(gradproj) < gradtol ) {
         keepgoing<-FALSE
         msg<-paste("Small gradient projection =",abs(gradproj),'')
         if (trace>0) cat(msg,"\n")
         break
      }
      ########################################################
      ####      Backtrack only Line search                ####
      changed<-TRUE # Need to set so loop will start
      steplength <- oldstep
      cat("steplength set to oldstep = ",steplength,"\n")
      while ((f >= fmin) && changed && (! accpoint)) {
         # We seek a lower point, but must change parameters too
         steplength<-min(steplength,maxstep)
         cat("steplength reset to ",steplength,"\n")
         bvec<-par+steplength*t
         if (trace>2) {
            cat("new bvec for steplength =",steplength,":")
            print(bvec)
         }
         changed <- ( ! identical((bvec+reltest), (par+reltest)) )
         if( changed ) { # compute new step, if possible
            f <- myfn(bvec, ...)
            ifn<-ifn+1
            if (ifn > maxfeval) {
               msg<-"Too many function evaluations"
               if (dowarn) warning(msg)
               termination<-1
               changed<-FALSE
               keepgoing<-FALSE
               break
            }
            if (attr(f,"inadmissible") ) {
               if(trace>0){
                 cat("Function is not calculable at bvec:")
                 print(bvec)
               }
               msg="Function is not calculable at an intermediate point"
	    }
            if (f < fmin) { # We have a lower point. Is it "low enough" i.e., acceptable
              accpoint <- ( f <= fmin + gradproj * steplength * acctol)
            } else {
              steplength<-steplength*stepredn
              if (trace > 0) cat("*")
            }
         } else { # NOT changed in step reduction
            if (trace>1) cat("Unchanged in step redn \n")
         } 
      }  # end while ((f >= fmin) && changed )
   } # end if gradproj<0
# === end backtrack - try qmin ===
## ?? If unchanged cannot have qmin best ??!!110506
   cat(msg,"\n")
   cat("f, fmin, gradproj, steplength: ",f, fmin, gradproj, steplength,"\n")
   if (f <= fmin + gradproj*steplength) { # cannot do qmin; f too low
      if (trace > 2) cat("f too low for quadratic min\n")
   } else {
      if (tryqmin) {
         f1<-f # save best fn
         step1<-steplength # to save these values
         if(trace>2) cat("f-backtrack OK =",f1," fmin=f0=",fmin," gradproj=",
                                      gradproj," step:",steplength,"\n")
         pstep <- 2.0*(f1 - fmin - gradproj * step1) 
         if (trace>2) cat("denominator =",pstep,"\n")
         if (pstep != 0.0) { pstep = -gradproj * step1 * step1 / pstep }
         if (trace > 2) cat("Raw pstep = ",pstep,"\n")
         pstep<-min(pstep,maxstep)
          if (pstep > 0) {
             bvec<-as.vector(par+pstep*t)
             changed <- ( ! identical((bvec+reltest), (par+reltest)) )
             if ( changed ) { # compute function, then pstep, if possible
                f2 <- myfn(bvec, ...) # This will fail if try() fails.
                ifn<-ifn+1
                if (trace > 0) cat("#")  # to indicate quadmin evaluation
                if (trace > 2) cat(" New fn value for pstep= ",pstep,"  f=",f2,"\n")
                if (f2<f1) { # update best
                   steplength<-pstep # no longer zero
                   f<-f2
                   msg<-"quadmin best"
                }
                if (ifn > maxfeval) { 
                   stop("Too many function evaluations")
		   ### ?? this NEEDS to be fixed so we exit properly
                }
             } else { # not changed
                msg<-"No change in parameters in quadratic fit stage of line search"
                if (trace > 0) cat(msg,"\n")
             } # end not changed
          } # end pstep>0
       } # end if tryqmin
    } # end quadfit not possible
# === end qmin ===
    if (accpoint) { # matrix update if acceptable point.
      if (bounds) {
        ## Reactivate constraints??
        for (i in 1:n) {
          if (bdmsk[i]==1) { # only interested in free parameters
            # make sure < not <= below to avoid Inf comparisons
            if ( (bvec[i]-lower[i]) < ceps*(abs(lower[i])+1.0) ) { # are we near or lower than lower bd
                if (trace > 2) cat("(re)activate lower bd ",i," at ",lower[i],"\n")
              bdmsk[i] <- -3
            } # end lower bd reactivate
            if ( (upper[i]-bvec[i]) < ceps*(abs(upper[i])+1.0) ) { # are we near or above upper bd
                if (trace > 2) cat("(re)activate upper bd ",i," at ",upper[i],"\n")
              bdmsk[i] <- -1
            } # end lower bd reactivate
          } # end test on free params
        } # end reactivate constraints
      } # if bounds
      g<-mygr(bvec, ...) # ?? use try()
      ig<-ig+1
      if (ig > maxit) {
         keepgoing=FALSE
         msg="Too many gradient evaluations"
         warning(msg)
         termination<-1
         break
      }
      t[which(bdmsk==0)] <- 0 # active mask or constraint (change to t not g 110506)
      t<-as.vector(steplength*t)
      c<-as.vector(g-c)
      D1<-sum(t*c)
      if (D1>0) { 
        y<-as.vector(B %*% c)
        D2<-as.double(1.0+(t(c) %*% y)/D1) # as.double because D2 is a 1 by 1 matrix otherwise
        # May be able to be more efficient below -- need to use outer function
        B <- B - (outer(t, y) + outer(y, t) - D2 * outer(t,t))/D1
      } else {
        if (trace>0) cat("UPDATE NOT POSSIBLE\n")
        ilast<-ig # note gradient evaluation when update failed
      } # D1 > 0 test 
      oldstep<-0.5*steplength+0.5
    } else { # no acceptable point
      if (trace > 0) cat("No acceptable point\n")
      if(ig == ilast) { # we reset to gradient and did new linesearch
        keepgoing<-FALSE # no progress possible
        if( termination < 0) {
           termination <- 0
           msg<-"Normal termination"
        }
        if(trace > 0) cat(msg,"\n")
      } else {
        ilast<-ig # reset to gradient search
        if (trace>0) cat("Reset to gradient search\n")
        oldstep<-1
      } # end else ig != ilast
    } # end else no accpoint 
  } # end main loop  (while keepgoing)
  if (trace > 0) cat("Seem to be done VM\n")
  ans<-list(par, fmin, c(ifn, ig), termination, msg, bdmsk)
## ?? need to fix up message
  names(ans)<-c("par", "value", "counts", "termination", "message", "bdmsk")
  return(ans)
} ## end of Rvmmin

