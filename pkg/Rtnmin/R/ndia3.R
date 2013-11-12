ndia3<-function(e, v, gv, r, vgv){
##---------------------------------------------------------
## update the preconditioning matrix based on a diagonal 
## version of the bfgs quasi-newton update.
##---------------------------------------------------------
  tol    <- 1e-6
  vr     <- as.numeric(crossprod(v,r))
#%   cat("ndia3: vgv, vr, then e:",vgv, vr,"\n")
#%   print(e)
  if (abs(vr)>tol && abs(vgv)>tol) {
    e <- e - (r*r)/vr + (gv*gv)/vgv # CAUTION! May be crossprod
    ind  <- which(e < tol)
    if (length(ind)>0){
       e[ind] <- 1
    }
  }
#%   print(e)
  e # Need to return the object?? Possibly not!
}