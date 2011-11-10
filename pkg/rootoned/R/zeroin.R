zeroin <- function(f, interval,
             tol = .Machine$double.eps^0.5, maxiter = 1000, trace=FALSE,...){
# trace <- FALSE # to allow for debugging
# epsilon <- .Machine$double.eps^0.5
epsilon <- .Machine$double.eps
if (trace) cat("epsilon= ",epsilon,"\n")

# uniroot <- function(f, interval, ...,
#             lower = min(interval), upper = max(interval),
#             f.lower = f(lower, ...), f.upper = f(upper, ...),
#             tol = .Machine$double.eps^0.25, maxiter = 1000)

#      a zero of the function  f(x)  is computed in the interval ax,bx .
#
#  input..
#
#  ax     left endpoint of initial interval
#  bx     right endpoint of initial interval
#  f      function subprogram which evaluates f(x) for any x in
#         the interval  ax,bx
#  tol    desired length of the interval of uncertainty of the
#         final result (.ge.0.)
#
#  output..
#
#  zeroin abscissa approximating a zero of  f  in the interval ax,bx

ax <- interval[1]
bx <- interval[2]
fa <- f(ax, ...)
fb <- f(bx, ...)
if (trace) cat("Start zeroin: f(",ax,")=",fa,"   f(",bx,")=",fb,"\n")

a <- ax
b <- bx
c <- a
fc <- fa;
wtol <- tol # to avoid changing tol
# First test if we have found a root at an endpoint 
    if(fa == 0.0) {
	wtol <- 0.0
	maxit <- 0
        if (trace) cat("fa is 0\n")
	return(list(root=a, froot=0, wtol=0.0, maxit=2))
    }
    if(fb ==  0.0) {
	wtol <- 0.0
	maxit <- 0
        if (trace) cat("fb is 0\n")
	return(list(root=b, froot=0, wtol=0.0, maxit=2))
    }
    maxit <- maxiter + 2 # count evaluations as maxiter-maxit

    while(maxit > 0)		## Main iteration loop	*/
    {
        if (trace) cat("Top of iteration, maxit=",maxit,"\n")

	prev_step <- b-a		## Distance from the last but one
					##   to the last approximation	*/
	#double tol_act;			## Actual tolerance		*/
	#double p;			## Interpolation step is calcu- */
	#double q;			## lated in the form p/q; divi-
	#				 * sion operations is delayed
	#				 * until the last moment	*/
	#double new_step;		## Step at this iteration	*/

	if( abs(fc) < abs(fb) )
	{				## Swap data for b to be the	*/
            if (trace) cat("fc smaller than fb\n")
	    a <- b
            b <- c
            c <- a	## best approximation		*/
	    fa<-fb
            fb<-fc
            fc<-fa
	}
	tol_act <- 2*epsilon*abs(b) + tol/2
        if (trace) cat("tol_act= ",tol_act,"\n")

	new_step <- (c-b)/2 # bisection
        if (trace) cat("new_step= ",new_step,"\n")

	if( (abs(new_step) <= tol_act) || (fb == 0.0 ))
	{
            if (trace) cat("DONE! -- small new_step or fb=0\n")
	    wtol = abs(c-b)
	    return(list(root=b, froot=fb, rtol=wtol, maxit=maxiter-maxit)) ## Acceptable approx. is found	*/
	}

	## Decide if the interpolation can be tried	*/
	if( (abs(prev_step) >= tol_act)	## If prev_step was large enough*/
	    && (abs(fa) > abs(fb)) ) {	## and was in true direction,
					## Interpolation may be tried	*/
	#    register double t1,cb,t2;
            if (trace) cat("prev_step larger than tol_act and fa bigger than fb\n")

	    cb <- c-b
	    if( a==c ) {		## If we have only two distinct	*/
					## points linear interpolation	*/
		t1 <- fb/fa		## can only be applied		*/
		p <- cb*t1
		q <- 1.0 - t1
                if (trace) cat("a == c\n")
	    }
	    else {			## Quadric inverse interpolation*/
                if (trace) cat("a != c\n")
		q <- fa/fc
                t1 <- fb/fc
                t2 <- fb/fa
		p <- t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) )
		q <- (q-1.0) * (t1-1.0) * (t2-1.0)
	    }
	    if( p>0.0 )	{	## p was calculated with the */
                if (trace) cat("p>0\n")
		q <- -q			## opposite sign; make p positive */
	    } else {			## and assign possible minus to	*/
                if (trace) cat("! p>0\n")
		p <- -p			## q				*/
            }
	    if( (p < (0.75*cb*q-abs(tol_act*q)/2))  ## If b+p/q falls in [b,c]*/
		&& (p < abs(prev_step*q/2)) ) {	## and isn't too large	*/
                if (trace) cat("p satisfies conditions for changing new_step\n")
		new_step <- p/q			## it is accepted
	    }					# * If p/q is too large then the
						# * bisection procedure can
						# * reduce [b,c] range to more
						# * extent */
	}

	if( abs(new_step) < tol_act) {	## Adjust the step to be not less*/
            if (trace) cat("new_step smaller than tol_act\n")
	    if( new_step > 0.0 )	## than tolerance		*/
		new_step <- tol_act
	    else
		new_step <- -tol_act
	}
	a <- b
        fa <- fb			## Save the previous approx. */
	b <- b + new_step
        fb <- f(b, ...)
        maxit <- maxit-1
	## Do step to a new approxim. */
	if( ((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0)) ) {
            if (trace) cat("make c to have sign opposite to b\n")
	    ## Adjust c for it to have a sign opposite to that of b */
	    c <- a
            fc <- fa
	}

    }
    ## failed! */
    if (trace) cat("Failed!\n")
    wtol <- abs(c-b)
    maxit <- -1
    return(list(root=b, froot=NA, rtol=wtol, maxit=maxit)) ## Acceptable approx. is found	*/
}

