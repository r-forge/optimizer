##################################################################
get.best <- function(optimx.obj, maximize=FALSE) {
    # returns solution that is best according to fvalue (i.e., minimum)
    # ?? for now ignore maximize
    # Let us take LAST value from optimx
    # optimx.obj = object returned by `optimx'
    if (maximize) fs<--1 else fs<-1
    fv<-unlist(optimx.obj$fvalues)
    ibest<-which(fv==min(fv*fs)) # fixed for maximize
    if (length(ibest)>1) ibest<-ibest[1] # in case of a tie
    if (length(ibest)<1) { # No solution
       ret<-list(par=NA, value=NA,counts=NA, conv=9999, message="No solution found",
             hessian=NA)
       return(ret)
    }
    par<-optimx.obj$par[[ibest]]
    value<-fv[ibest]
    counts<-c(optimx.obj$fns[[ibest]], optimx.obj$grs[[ibest]])
    convergence<-optimx.obj$conv[[ibest]]
    message<-attr(optimx.obj, "details")[[ibest]]$message
    hessian<-attr(optimx.obj, "details")[[ibest]]$nhatend
    ret<-list(par=par,
              value=value,
              counts=counts,
              convergence=convergence,
              message=message,
              hessian=hessian)
}
##################################################################
