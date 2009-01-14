
BBsolve <- function(par, fn, algorithm = "dfsane", methods = c(2,1), 
     M = c(10, 100), NM = c(FALSE, TRUE) ) {
control.pars <- expand.grid(methods, M, NM)
ans.best.value <- Inf
ans <- vector("list", length=nrow(control.pars))

for (i in 1: nrow(control.pars) ) {
   cpars <- unlist(control.pars[i, ])
   cat("Try : ", i, "method = ", cpars[1],  "M = ", cpars[2],   "NM  = ", as.logical(cpars[3]), "\n")
   if (algorithm == "dfsane")
       ans[[i]] <- dfsane(par, fn, method=as.numeric(cpars[1]), 
                      control=list(M=as.numeric(cpars[2]), NM=as.logical(cpars[3]), trace=FALSE))
   if (algorithm == "sane")
       ans[[i]] <- sane(par, fn, method=as.numeric(cpars[1]), 
           control=list(M=as.numeric(cpars[2]), NM=as.logical(cpars[3]), trace=FALSE))

   if (ans[[i]]$residual < ans.best.value) {
       ans.best <- ans[[i]]
       ans.best.value <- ans.best$residual
       }
   if (ans[[i]]$convergence  == 0) break
   }
if (ans[[i]]$convergence != 0)
    warning("Unsuccessful convergence. Try again with `SANE' or repeat with a different starting value.\n Here is the best solution:\n")
ans.best
}
