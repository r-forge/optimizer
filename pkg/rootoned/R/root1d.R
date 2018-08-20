root1d <- function(f, interval,
                   tol = .Machine$double.eps^0.5, maxiter = 1000, trace=FALSE, ...) {
  
  # trace <- FALSE
  #     uniroot(f, interval,
  #             lower = min(interval), upper = max(interval),
  #             f.lower = f(lower, ...), f.upper = f(upper, ...),
  #             tol = .Machine$double.eps^0.25, maxiter = 1000, ...)
  
  #var lbound, ubound: real;
  #                 var ifn: integer;
  #                     tol : real;
  #                 var noroot: boolean );
  
  #var
  # nbis: integer;
  # b, fb, flow, fup : real;
  # notcomp: boolean;
  
  
  if (trace) cat('alg18 == root1d -- root of a function of one variable\n')
  
  ubound <- interval[2]
  lbound <- interval[1]
  
  #  notcomp := FALSE; -- use try()
  ifn <- 2
  nbis <- 5
  fup <- f(ubound,...)
  fup0<-fup # save it
  #  if notcomp then halt;
  flow <- f(lbound,...)
  flow0<-flow # save it
  #  if notcomp then halt;
  #if (trace) cat('f(',lbound,')=',flow,'  f(',ubound,')=',fup,"\n")
  
  if (sign(fup)*sign(flow) > 0) {
    noroot  <-  TRUE
  } else { noroot  <-  FALSE }
  op<-"start"
  lastwidth <- (ubound - lbound)
  while ( (! noroot) && ((ubound-lbound)>tol) && (ifn < maxiter)) {
    if (trace) cat("f(",lbound,")=",flow,"  f(",ubound,")=",fup,"  interval=",ubound-lbound,"  ",op,"\n")
    cat("Width reduced: ", (ubound - lbound  >= lastwidth - tol))
    if ( (( (nbis * ((ifn - 2) %/% nbis) == (ifn - 2)) ) && (ifn > 2)) || (ubound - lbound  >= lastwidth - tol) ){
      op<-"Bisect"
      #if (trace) cat('Bisect  \n')
      b  <-  lbound + 0.5*(ubound - lbound)
    } else {
      #if (trace) cat('False position \n')
      op<-"FalsePos"
      ## 180820 JN: May give trouble if magnitudes wild
      b  <-  (lbound*fup-ubound*flow)/(fup-flow)
    }
    
    #  cat("b =",b,"  lbound=",lbound,"  ubound=",ubound,"\n")
    ## Edited by Oliver Dechant July 30th, 2018
    ## 180820 Temporarily comment out. Do these stmts work for scalars?
    #  b[is.na(b)] <- 0
    #  lbound[is.na(lbound)] <- 0
    ## END of edits
    
    #    if (b <= lbound) {
    #      b  <-  lbound
    #      ubound  <-  lbound # force convergence? We cannot compute an in-range value.
    #    }
    #    if (b >= ubound) { 
    #      b  <-  ubound
    #      lbound  <-  ubound # force convergence? We cannot compute an in-range value.
    #    } # But we DO NOT break successfully -- not following code in CNM79
    if (b <= lbound) {
      op <- "Fail-FP-low"
      if (b <= lbound) {
        ubound <- lbound
      } else { lbound <- ubound} # close interval
      cat("Fail-FP\n")
      break
    }
    ifn  <-  ifn+1
    fb  <-  f(b, ...)
    #    if notcomp then halt; ## not implemented yet (would need try())
    if (trace) cat(ifn,' evalns: f(',b,')=',fb,"\n")
    if (trace) cat('  width interval= ',(ubound-lbound),"\n")
    if ( (ubound-lbound)>tol ){
      if (sign(fb)*sign(flow) < 0.0) {
        fup  <-  fb
        ubound  <-  b
      } else {
        flow  <-  fb
        lbound  <-  b
      }
      lastwidth <- ubound - lbound
    }
  } # end while
  if (noroot) { 
    b<-NA # We did not find a root
    fb<-NA
  } # end while loop
  if (trace)   cat('Terminated at f(',b,')=',fb,"\n")
  if (trace)   cat('  Final interval width =',ubound-lbound,"\n")
  if (! is.na(fb) && (abs(fb) > 0.5*max(abs(fup0), abs(flow0)))) {
    msg<-"Final function magnitude > 0.5 * max(abs(f(lbound)), abs(f(ubound)))"
    if (trace) cat(msg,"\n")
    warning(msg)
  }
  res<-list(root=b, froot=fb, rtol=(ubound-lbound), fcount=ifn)
}


