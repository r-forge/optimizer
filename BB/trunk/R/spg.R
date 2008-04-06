spg <- function(par, fn, gr=NULL, method=3, project=NULL, lower=-Inf, upper=Inf,  control=list(),  ... ) {
############################################################
# Non-monotone spectral projected-gradient method for minimization
# Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected gradient methods on convex sets, 
# SIAM J Optimization, 10, 1196-1211.
# Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for convex-constrained optimization, 
# ACM Transactions on Mathematical Software.
###############################################
# R adaptation, with significant modifications, by  Ravi Varadhan, Johns Hopkins University, March 25, 2008.
#
#   Most important modification is the availability of different options for Barzilai-Borwein steplengths
#   Three different Barzilai-Borwein steplength options can be chosen.
#   Method = 1 is the steplength used in Birgin EG, Martinez JM, and Raydan M (2000)  
#   Method = 2 is another BB steplength proposed in Barzilai and Borwein's (1988) original paper 
#   Method = 3, is a new steplength, first proposed in Varadhan and Roland (2008).
#
#   Method = 3 is the "default" since it performed slightly better than others in our numerical experiments in terms of convergence to better optimum.
#
# Please refer to Varadhan and Gilbert (2008, unpublished) for details
#
# Incorporates box constraints
# The user can define his/her own projection function "project" to handle more complicated constraints
#
################################################
    ctrl <- list(M=10, maxit=1500, ftol=1e-08, gtol=1.e-04, maxfeval=10000, maximize=FALSE, trace=TRUE, 
    triter=10, grad.method="simple", eps=1.e-07) # defaults
    ctrl[names(control)] <- control
    M     <- ctrl$M
    maxit <- ctrl$maxit
    ftol  <- ctrl$ftol
    gtol  <- ctrl$gtol
    maxfeval <- ctrl$maxfeval
    maximize <- ctrl$maximize
    trace <- ctrl$trace
    triter <- ctrl$triter
    grad.method <- ctrl$grad.method
    eps <- ctrl$eps
###########################################################
nmls <- function(p, f, d, gtd, lastfv, feval, func, maxfeval, ... ){
#  local function
# Non-monotone line search of Grippo with safe-guarded quadratic interpolation
gamma <- 1.e-04
fmax <- max(lastfv)
alpha <- 1
pnew <- p + alpha*d
fnew <- try(func(pnew , ...),silent=TRUE)
feval <- feval + 1

if (class(fnew)=="try-error" | is.nan(fnew)) return(list(p=NA, f=NA, feval=NA, lsflag=1))

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

    if (class(fnew)=="try-error" | is.nan(fnew)) return(list(p=NA, f=NA, feval=NA, lsflag=1))


    if (feval > maxfeval) return(list(p=NA, f=NA, feval=NA, lsflag=2))

}  # Main while condition loop ends

return(list(p=pnew, f=fnew, feval=feval, lsflag=0))
}

#############################################
if (is.null(project)) {
project.box <- function(x, lower, upper) {
# local function
# Defined only when user-doesn't specify his/her own projection algorithm
# Projecting to ensure that box-constraints are satisfied
#

x[x < lower] <- lower[x < lower]
x[x > upper] <- upper[x > upper]
return(x)
}
}
#############################################

#  Initialization
lmin <- 1.e-30
lmax <- 1.e30
iter <- 0
feval <- 0
geval <- 0
lastfv <- rep(-1.e99, M)
fbest <- NA
fchange <- fchg.rel <- 1

if (maximize) func <- function(x, ...) {-1 * fn(x, ...)} 
else func <- function(x, ...) fn(x, ...)

# Project initial guess
if (is.null(project)) par <- try(project.box(par, lower, upper), silent=TRUE) 
else par <- try(project(par, ...), silent=TRUE)

if (class(par) == "try-error" | any(is.nan(par)) | any(is.na(par))) {
cat("\n Failure: Error in projecting initial guess! \n ")
return(NULL)
} else pbest <- par

f <- try(func(par, ...),silent=TRUE)      
feval <- feval + 1
if (class(f)=="try-error" | is.nan(f) | is.na(f)){
cat("\n Failure: Error in initial function evaluation! \n")
return(NULL)
}
f0 <- fbest <- f


if (is.null(gr)) require(numDeriv)

if (is.null(gr)) g <- try(grad(func=func, x=par, method=grad.method, method.args=list(eps=eps, r=2), ...),silent=TRUE)
else g <- try(gr(par, ...),silent=TRUE)
geval <- geval + 1


if (class(g)=="try-error" | any(is.nan(g)) | any(is.na(g)) ){
cat("\n Failure: Error in initial gradient evaluation! \n")
return(NULL)
}

lastfv[1] <- f
fbest <- f
pg <- par - g

if (is.null(project)) pg <- project.box(pg, lower, upper)
else pg <- project(pg, ...)
pg <- pg - par

pg2n <- sqrt(sum(pg*pg))
pginfn <- max(abs(pg))
gbest <- pg2n
if (pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))

continue <- TRUE

#######################
#  Main iterative loop
#######################
while(continue) {
iter <- iter + 1
d <- par - lambda * g

if (is.null(project)) d <- try(project.box(d, lower, upper), silent=TRUE)
else d <- try(project(d, ...), silent=TRUE)

if (class(d) == "try-error" | any(is.nan(d)) | any(is.na(d))) {
lsflag <- 4
break
} 

d <- d - par
gtd <- sum(g * d)

nmls.ans <- nmls(par, f, d, gtd, lastfv, feval , func, maxfeval, ...)
lsflag <- nmls.ans$lsflag

if(lsflag != 0) break
 
fchange <- abs(f - nmls.ans$f)
fchg.rel <- fchange / abs(f)
f <- nmls.ans$f 
pnew <- nmls.ans$p
feval <- nmls.ans$feval
lastfv[(iter %% M) + 1] <- f

if (is.null(gr)) gnew <- try(grad(func=func, x=pnew, method=grad.method, method.args=list(eps=eps, r=2), ...),silent=TRUE)
else gnew <- try(gr(pnew, ...),silent=TRUE)
geval <- geval + 1

if (class(gnew)=="try-error" | any(is.nan(gnew))){
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
g <- gnew
pg <- par - g

if (is.null(project)) pg <- try(project.box(pg, lower, upper), silent=TRUE)
else pg <- try(project(pg, ...), silent=TRUE)

if (class(pg) == "try-error" | any(is.nan(pg)) | any(is.na(pg))) {
lsflag <- 4
break
} 

pg <- pg - par
pg2n <- sqrt(sum(pg*pg))
pginfn <- max(abs(pg))
ginfn <- max(abs(g))

f.rep <- (-1)^maximize * f
if (trace & (iter%%triter == 0)) cat("iter: ",iter, " f-value: ", f.rep, " grad: ",pginfn, "\n")

if (f < fbest) {
fbest <- f
pbest <- pnew
gbest <- pginfn
}

continue <- (fchange > ftol) & (fchg.rel > ftol) & (pginfn > gtol) & (iter <= maxit)

}   # while condition loop concludes

if (lsflag==0) {
if (fchange <= ftol | fchg.rel <= ftol) conv <- list(type=0, message="Successful convergence")
if (iter >= maxit) conv <- list(type=1, message="Maximum number of iterations exceeded")
if (pginfn <= gtol & ginfn > gtol) conv <- list(type=3, message="Only projected gradient convergence")
} else {
	par <- pbest
	f.rep <- f <- (-1)^maximize * fbest
	pginfn <- gbest
if (lsflag==1) conv <- list(type=4, message="Failure:  Error in function evaluation")
if (lsflag==2) conv <- list(type=2, message="Maximum function evals exceeded")
if (lsflag==3) conv <- list(type=5, message="Failure:  Error in gradient evaluation")
if (lsflag==4) conv <- list(type=6, message="Failure:  Error in projection")
}

absred <- (-1)^maximize * (f0 - f)
relred <- 1 - ((f0-f)/f0)
ginfn <- max(abs(g))
browser()
return(list(par=par, value=f.rep, abs.reduction=absred, rel.reduction=relred, grad=pginfn, iter=iter, feval=feval,  
    convergence=conv$type, message=conv$message))
}

