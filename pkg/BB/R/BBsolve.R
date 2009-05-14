BBsolve <- function(par, fn, algorithm = c("dfsane", "sane"), 
	method=c(2,1,3), M = c(10, 100), NM=c(TRUE, FALSE), ... , 
	tol=1.e-07, maxit=1500) 
{
if (is.null(dim(par))) par <- matrix(par, nrow=1)
if (ncol(par) > 20) NM <- FALSE   

control.pars <- expand.grid(method=method, M=M, NM=NM)
ans.best.value <- Inf
ans <- vector("list", length=nrow(par))
alg <- match.arg(algorithm)

feval <- iter <- s <- 0
pmat <- matrix(NA, nrow(par) ,ncol(par))

for (k in 1:nrow(par)){
cat("Parameter set : ", k, "\n")
	p <- par[k, ]

for (i in 1: nrow(control.pars) ) {
   cpars <- unlist(control.pars[i, ])
   #cat("Try : ", i, "Method = ", cpars[1], "M = ", cpars[2], "Nelder-Mead = ", cpars[3], "\n")

   if (alg == "dfsane") temp <- 
     dfsane(par=p, fn, method=cpars[1], control=list(M=as.numeric(cpars[2]), 
        NM=cpars[3], 
	maxit=maxit, tol=tol, trace=FALSE, noimp=min(100, 5*cpars[2])), ...)

   if (alg == "sane") temp <-
       sane(par=p, fn, method=cpars[1], control=list(M=as.numeric(cpars[2]),
        NM=cpars[3], 
	maxit=maxit, tol=tol, trace=FALSE, noimp=min(100, 5*cpars[2])), ...)

   feval <- feval + temp$feval
   iter <- iter + temp$iter

   if (temp$convergence  == 0) {
	cat ("Successful convergence \n\n")
	s <- s + 1
	ans[[k]] <- temp
	pmat[k, ] <- temp$par
	ans[[k]]$cpar <- cpars
	break
	} else if (temp$residual < ans.best.value) {
	ans.best <- temp
	ans.best.value <- ans.best$residual
	ans.best$feval <- feval
	ans.best$iter <- iter
	ans.best$cpar <- cpars
	}

}  # "i" loop completed

if (temp$convergence != 0) cat ("Unsuccessful convergence \n\n")

}  # "k" loop completed

if (s == 0) ans[[1]] <- ans.best
success <- !is.na(pmat[,1])

if (k ==1 | s == 0) ans[[1]] else list(par = pmat[success, ], info=ans[success])
}

