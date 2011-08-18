root1d <- function(f, interval,
             tol = .Machine$double.eps^0.5, maxiter = 1000, ...) {

debug <- FALSE
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


if (debug) cat('alg18 == root1d -- root of a function of one variable\n')

ubound <- interval[2]
lbound <- interval[1]

#  notcomp := FALSE; -- use try()
  ifn <- 2
  nbis <- 5
  fup <- f(ubound,...)
#  if notcomp then halt;
  flow <- f(lbound,...)
#  if notcomp then halt;
if (debug) cat('f(',lbound,')=',flow,'  f(',ubound,')=',fup,"\n")

  if (fup*flow>0) {
     noroot  <-  TRUE
  } else { noroot  <-  FALSE }

  while ( (! noroot) && ((ubound-lbound)>tol) ) {

    if ( (nbis * ((ifn - 2) %/% nbis) == (ifn - 2)) ){
if (debug) cat('Bisect  \n')
      b  <-  lbound + 0.5*(ubound - lbound)
    } else {
if (debug) cat('False position \n')
      b  <-  (lbound*fup-ubound*flow)/(fup-flow)
    }

    if (b <= lbound) {
      b  <-  lbound
      ubound  <-  lbound
    }
    if (b >= ubound) { 
      b  <-  ubound
      lbound  <-  ubound
    }
    ifn  <-  ifn+1

    fb  <-  f(b, ...)
#    if notcomp then halt;
if (debug) cat(ifn,' evalns: f(',b,')=',fb,"\n")
if (debug) cat('  width interval= ',(ubound-lbound),"\n")
    if ( (ubound-lbound)>tol ){
      if (fb*flow<0.0) {
        fup  <-  fb
        ubound  <-  b
      } else {
        flow  <-  fb
        lbound  <-  b
      }
    }
  } # end while
if (debug)   cat('Converged to f(',b,')=',fb,"\n")
if (debug)   cat('  Final interval width =',ubound-lbound,"\n")
  res<-list(root=b, froot=fb, rtol=(ubound-lbound), fcount=ifn)
}
