# nlmrt.R -- print and summary methods
#   result <- list(resid = resbest, jacobian = Jac, feval = feval, 
#        jeval = jeval, coeffs = pnum, ssquares = ssbest)
#    class(result) <- "nlmrt"

print.nlmrt <- function(x, ...) {
    smalltol <- .Machine$double.eps * 1000
    options(digits = 7) # ??this is default
    resname <- deparse(substitute(x))
    cat("nlmrt class object:",resname,"\n")
    coef <- x$coeffs
    pname<-names(coef)
    npar <- length(coef)
    lo <- x$lower
    if (is.null(lo)) lo <- rep( -Inf, npar)
    up <- x$upper
    if (is.null(up)) up <- rep( Inf, npar)
    mi <- x$maskidx
    mt <- rep(" ",npar) # start with all "unmasked"
    mt[mi] <- "M" # Put in the masks
    ct <- rep(" ",npar) # start with all "free"
    for (i in seq_along(coef)){
       if (lo[[i]] - coef[[i]] > 0) {
          ct[[i]] <- "-" # lower bound violation
       } else { 
          if (coef[[i]] - lo[[i]] < smalltol*(abs(coef[[i]])+smalltol) ) ct[[i]] <- "L" # "at" lower bound
       }
       if (coef[[i]] - up[[i]] > 0) {
          ct[[i]] <- "+" # lower bound violation
       } else { 
          if (up[[i]] - coef[[i]] < smalltol*(abs(coef[[i]])+smalltol) ) ct[[i]] <- "U" # "at" upper bound
       }
    }
#    cat("  name    ","  coeff      ","    SE    "," gradient ","  tstat  ","  pval   ","\n")
    cat("  name    ","  coeff      ","\n")
    # separate entries so they can be played with      
#    cat("coeffs:")
    for (i in seq_along(coef)){
       pc <- pname[[i]]
       # ?? adjust length of pc to 8 ?? chars (not needed in print, but in summary
       # add index i??
       cat(format(pc, width=8),"  ")
       cat(format(coef[[i]], width=8)," ")
       cat(ct[[i]],mt[[i]]," ")
       cat("\n")
    }
    cat("ssquares = ",x$ssquares,"\n")
    invisible(x)
}
# working 120901
## add
#  bounds and masks

summary.nlmrt <- function(object, ...) {
    smalltol <- .Machine$double.eps * 1000
    options(digits = 5) # 7 is default
    resname <- deparse(substitute(object))
    cat("Summary of nlmrt class object ",resname,"\n")
#    cat("coeffs:")
#    print(object$coeffs)
#    cat("ssquares = ",object$ssquares,"\n")
    JJ <- object$jacobian
    res <- object$resid
    coef <- object$coeffs
    resname <- deparse(substitute(x))
    cat("nlmrt class object:",resname,"\n")
    pname<-names(coef)
    npar <- length(coef)
    lo <- object$lower
    if (is.null(lo)) lo <- rep( -Inf, npar)
    up <- object$upper
    if (is.null(up)) up <- rep( Inf, npar)
    mi <- object$maskidx
    mt <- rep(" ",npar) # start with all "unmasked"
    mt[mi] <- "M" # Put in the masks
    bdmsk <- rep(1, npar) # bounds and masks indicator ?? should it be 1L
    bdmsk[mi] <- 0 # masked
    cat("bdmsk:")
    print(bdmsk)
    ct <- rep(" ",npar) # start with all "free"
    for (i in seq_along(coef)){
       if (lo[[i]] - coef[[i]] > 0) {
          ct[[i]] <- "-" # lower bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -3
       } else { 
          if (coef[[i]] - lo[[i]] < smalltol*(abs(coef[[i]])+smalltol) ) {
             ct[[i]] <- "L" # "at" lower bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -3 # leave mask indication intact
          }
       }
       if (coef[[i]] - up[[i]] > 0) {
          ct[[i]] <- "+" # lower bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -1
       } else { 
          if (up[[i]] - coef[[i]] < smalltol*(abs(coef[[i]])+smalltol) ) {
             ct[[i]] <- "U" # "at" upper bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -1 # leave mask indication intact
          }
       }
    }
    cat("ct:")
    print(ct)
    cat("new bdmsk:")
    print(bdmsk)
    ss <- object$ssquares
    nobs <- length(res)
    ndof <- nobs - npar
    if (ndof <= 0) stop(paste("Inadmissible degrees of freedom =",ndof,sep=''))
    sighat2 <- ss/(ndof)
    dec <- svd(JJ)
    U <- dec$u
#    cat("U:")
#    print(U)
    V <- dec$v
#    cat("V:")
#    print(V)
    Sd <- dec$d
#    cat("Sd:")
#    print(Sd)
#    if (min(Sd) <= smalltol * max(Sd)) { # singular
#       SEs <- rep(Inf, npar)
#    } else {
       Sinv <- 1/Sd
       Sinv[which(bdmsk != 1)] <- 0
       cat("Sinv:")
       print(Sinv)
       VS <- crossprod(t(V), diag(Sinv))
       cat("VS:")
       print(VS)
       Jinv <- crossprod(t(VS))
       var <- Jinv * sighat2
       SEs <- sqrt(diag(var))
#    }
    gr <- crossprod(JJ, res)
    tstat <- coef/SEs
    pval<-2*(1-pt(tstat, df=ndof))
    # object digits??
    cat("  name    ","  coeff        ","    SE    "," gradient   ","  tstat  ","   pval   ","\n")
    for (i in seq_along(coef)){
       pc <- pname[[i]]
       # ?? adjust length of pc to 8 ?? chars (not needed in print, but in summary
       # add index i??
       cat(format(pc, width=8),"  ")
       cat(format(coef[[i]], width=8)," ")
       cat(ct[[i]],mt[[i]],"  ")
       cat(format(SEs[[i]], width=8)," ")
       cat(format(gr[[i]], width=11)," ")
       cat(format(tstat[[i]], width=8)," ")
       cat(format(pval[[i]], digits=4, width=8)," ")
       cat("\n")
    }
    invisible(object)
}

## To add 
#    gradients
#    tstats
#    significance
#    bounds/masks -- will require return in nlmrt class objects

