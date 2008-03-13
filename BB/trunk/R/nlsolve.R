##########################################
nlsolve <- function(par, fn, lower=-Inf, upper=Inf, control=list(), ...) {
###########################
#  A solver for finding a zero of a system of non-linear equations
#  Actually minimizes the squared-norm of the set of functions by calling optim()
#  Uses the L-BFGS-B algorithm within optim()
#  All the control parameters can be passed as in the call to optim()
#
#  Author:  Ravi Varadhan, Center on Aging and Health, Johns Hopkins University, #rvaradhan@jhmi.edu
#  June 21, 2007
    ctrl <- list(nstarts=10) # defaults
    ctrl[names(control)] <- control
    nstarts     <- ctrl$nstarts

func <- function(x, ...) sum(fn(x, ...)^2)

ans <- try(optim(par, fn=func, method="L-BFGS-B", lower=lower, upper=upper,
...), silent=TRUE)

if (class(ans) != "try-error") err <- sqrt(ans$val/length(par))  
else err <- Inf

p0 <- par
istart <- 1
ntry <- 1

# Trying more random starts if the root was not located in the first try
while (err > 0.01 & istart < nstarts & ntry < 10*nstarts) {
par <- p0 * ( 1 + runif(length(par), -1, 1) )  # generating random starting values
#  Ensuring that constraints are satisfied
par[par < lower] <- lower[par < lower]
par[par > upper] <- upper[par > upper]

ans <- try(optim(par, fn=func, method="L-BFGS-B", lower=lower, upper=upper, ...), silent=TRUE)
ntry <- ntry + 1
if (class(ans) != "try-error") {
err <- sqrt(ans$val/length(par))  
istart <- istart + 1
}
}

if (err > 0.01) cat(" \n Caution:  Zero of the function has not been located!!! \n Increase nstart \n \n")

ans$value <- sqrt(ans$value / length(par))
ans
}

#
##########################################
