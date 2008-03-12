spg <- function(p, func, grad=NULL, lower=-Inf, upper=Inf, M=10, ftol=1.e-08, gtol=1.e-04, maxit=2500, maxfeval=10000, 
	trace=TRUE , grad.method="simple", ...) {
############################################################
# Non-monotone spectral projected-gradient method for minimization
# Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected gradient methods on convex sets, 
# SIAM J Optimization, 10, 1196-1211.
# Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for convex-constrained optimization, 
# ACM Transactions on Mathematical Software.
###############################################
# R translation (with minor modifications):  Ravi Varadhan, Johns Hopkins University, February 23, 2008.
# Incorporates box constraints
################################################
 
###########################################################
#  local function
nmls <- function(p, f, d, gtd, lastfv, feval, func, maxfeval, ... ){
# Non-monotone line search of Grippo with safe-guarded quadratic interpolation
gamma <- 1.e-04
fmax <- max(lastfv)
alpha <- 1
pnew <- p + alpha*d
fnew <- try(func(pnew , ...),silent=TRUE)
feval <- feval + 1

while (class(fnew)=="try-error" | is.nan(fnew)){
	if (alpha <= min(1.e-07, 1.e-04/(1 + sqrt(sum(d*d)))) ) return(NULL)

alpha <- alpha/4
pt <- p + alpha*d
fnew <- try(func(pt, ...),silent=TRUE)
feval <- feval + 1
}

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

	if (class(fnew)=="try-error" | is.nan(fnew)) return(NULL)
	if (feval > maxfeval) return(NULL)

}  # Main while condition loop ends

return(list(p=pnew, f=fnew, feval=feval))
}

#############################################
# local function
project <- function(x, lower, upper) {
# Projecting to ensure that box-constraints are satisfied

x[x < lower] <- lower[x < lower]
x[x > upper] <- upper[x > upper]
return(x)
}
#############################################


#  Initialization
lmin <- 1.e-30
lmax <- 1.e30
iter <- 0
feval <- 0
geval <- 0
lastfv <- rep(-1.e99, M)
pbest <- p
fbest <- NA
stagn <- FALSE
nsy <- 0
fchange <- 1

if (is.null(grad)) require(numDeriv)

f <- try(func(p, ...),silent=TRUE)      
feval <- feval + 1
if (class(f)=="try-error" | is.nan(f)){
cat("\n Failure: Error in function evaluation! \n Try another starting value \n")
return(NULL)
}

if (is.null(grad)) g <- try(grad(func=func, p, method=grad.method, method.args=list(r=2), ...),silent=TRUE)
else g <- try(grad(p, ...),silent=TRUE)
geval <- geval + 1

if (class(g)=="try-error" | any(is.nan(g))){
cat("\n Failure: Error in gradient evaluation! \n Try another starting value \n")
return(NULL)
}

lastfv[1] <- f
fbest <- f
pg <- p - g
pg <- project(pg, lower, upper)
pg <- pg - p

pg2n <- sqrt(sum(pg*pg))
pginfn <- max(abs(pg))
gbest <- pg2n
if (pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))

#######################
#  Main iterative loop
#######################
while( (fchange > ftol | pg2n > gtol) & iter <= maxit & !stagn) {
iter <- iter + 1

d <- p - lambda * g
d <- project(d, lower, upper)

d <- d - p
gtd <- sum(g * d)

nmls.ans <- nmls(p, f, d, gtd, lastfv, feval , func, maxfeval, ...)
if(is.null(nmls.ans)) {
cat("\n Failure: Error in line-search! \n Try another starting value \n")
return(NULL)
}
 
fchange <- abs(f - nmls.ans$f)
f <- nmls.ans$f 
pnew <- nmls.ans$p
feval <- nmls.ans$feval
lastfv[(iter %% M) + 1] <- f

if (is.null(grad)) gnew <- try(grad(func=func, pnew, method=grad.method, method.args=list(r=2), ...),silent=TRUE)
else gnew <- try(grad(pnew, ...),silent=TRUE)
geval <- geval + 1

if (class(gnew)=="try-error" | any(is.nan(gnew))){
cat("\n Failure: Error in gradient evaluation! \n Try another starting value \n")
return(NULL)
}

s <- pnew - p
y <- gnew - g
sts <- sum(s*s)
yty <- sum(y*y)
sty <- sum(s*y)

if (abs(sty) < 1.e-20) stagn <- TRUE
lambda <- min(lmax, max(lmin, sts/abs(sty)))

p <- pnew
g <- gnew
pg <- p - g
pg <- project(pg, lower, upper)

pg <- pg - p
pg2n <- sqrt(sum(pg*pg))
pginfn <- max(abs(pg))

if (trace & iter%%10 == 0) cat("iter: ",iter, " func: ", f, " grad: ",pg2n, "\n")

if (f < fbest) {
	fbest <- f
	pbest <- pnew
	gbest <- pg2n
	}

}   # while condition loop concludes

if (fchange <= ftol | pg2n <= gtol) conv <- list(type=0, message="Successful convergence")
if (stagn) conv <- list(type=1, message="Partial convergence due to Stagnation")
if (iter >= maxit) conv <- list(type=2, message="Maximum number of iterations exceeded")
 
#if (is.null(grad)) detach("package:numDeriv") 
return(list(par=pbest, func=fbest, grad=gbest, iter=iter, feval=feval, geval=geval, 
	convergence=conv$type, message=conv$message))
}


