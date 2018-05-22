### -*- mode: R; eval: (ess-set-style 'RStudio) -*-
##############################################
## BugReport_nlm_SourceCode.R
##  by Marie Boehnstedt
##    March 31, 2017
##  added 'neworder' argument: Martin Maechler
##############################################
##
###############################################
##
## Example 1: Rosenbrock banana valley function
##
f.Rosenb <- function(x1, x2) 100*(x2 - x1*x1)^2 + (1-x1)^2
grRosenb <- function(x1, x2) c(-400*x1*(x2 - x1*x1) - 2*(1-x1), 200*(x2 - x1*x1))
hessRosenb <- function(x1, x2) {
  a11 <- 2 - 400*x2 + 1200*x1*x1
  a21 <- -400*x1
  matrix(c(a11, a21, a21, 200), 2, 2)
}

fg <- function(x) {   # analytic gradient only
  x1 <- x[1]; x2 <- x[2]
  structure(f.Rosenb(x1, x2), "gradient" = grRosenb(x1, x2))
}
##
fgh <- function(x) {  # analytic gradient and Hessian
  x1 <- x[1];  x2 <- x[2]
  structure(f.Rosenb(x1, x2),
            "gradient" = grRosenb(x1, x2),
            "hessian" = hessRosenb(x1, x2))
}

system.time(for(i in 1:10) nlm.fg  <- nlm(fg,  c(-1.2, 1)) ) # 0.007
system.time(for(i in 1:10) nlm.fgh <- nlm(fgh, c(-1.2, 1)) ) # 0.009
str(nlm.fg)  ## ok: solution = [1 1]
str(nlm.fgh) ## wrong, but ok once  MM has fixed choldc() C code inside nlm() !
## but even then: both taking 24 iterations, and 'fg' is even slightly more accurate:
rbind(fg = nlm.fg $estimate,
      fgh= nlm.fgh$estimate) - 1
## fg   9.063217e-11  1.752587e-10
## fgh -7.145425e-10 -1.757105e-09
#############
##
##
## Example 2: Wood function
##
##' Wood function (length-4 vector 'x')
f.wood <- function(x) {
  stopifnot(is.numeric(x), length(x) == 4)
  100*(x[1]^2-x[2])^2 + (1-x[1])^2 + 90*(x[3]^2-x[4])^2 + (1-x[3])^2 +
    10.1*((1-x[2])^2 + (1-x[4])^2) + 19.8*(1-x[2])*(1-x[4])
}
##' Wood function (4 arguments 'x1' ... 'x4')
fwood <- function(x1,x2,x3,x4) {
  100*(x1^2-x2)^2 + (1-x1)^2 + 90*(x3^2-x4)^2 + (1-x3)^2 +
    10.1*((1-x2)^2 + (1-x4)^2) + 19.8*(1-x2)*(1-x4)
}
## automatically construct correct gradient and hessian:
woodf.gh <- function(x) {
  stopifnot(is.numeric(x))
  woodGH <- deriv3(body(fwood)[[2]],
                   c("x1","x2","x3","x4"), function.arg=TRUE)
  if(length(x) == 4)
    woodGH(x[1],x[2],x[3],x[4])
  else if(is.matrix(x) && ncol(x) == 4)
    woodGH(x[,1], x[,2], x[,3], x[,4])
  else stop("'x' must have length 4 or be a matrix with 4 columns")
}


##' gradient [wood function]
wfg <- function(x) {
  stopifnot(is.numeric(x), length(x) == 4)
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  c(g1,g2,g3,g4)
}

##' hessian [wood function]
wfh <- function(x) {
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2;     h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  matrix(c(h11,h12,h13,h14,
           h12,h22,h23,h24,
           h13,h23,h33,h34,
           h14,h24,h34,h44),ncol = 4)

}

## wood function with analytic gradient only:
woodfunc.g <- function(x) structure(f.wood(x), gradient = wfg(x))
#### wood function with analytic gradient and Hessian:
woodfunc.gh <- function(x) structure(f.wood(x), gradient = wfg(x), hessian = wfh(x))

nlm.wfg  <- nlm(p = c(-3,-1,-3,-1), woodfunc.g)
nlm.wfgh <- nlm(p = c(-3,-1,-3,-1), woodfunc.gh)
nlm.wfgh2<- nlm(p = c(-3,-1,-3,-1), woodf.gh) # using the automatically constructed {g,h}
str(nlm.wfg) # fine
str(nlm.wfgh)# not ok
str(nlm.wfgh2)# not ok
all.equal(nlm.wfgh, nlm.wfgh2)# they are the same (and both wrong "the same")
## "==>" woodf.gh() and woodfunc.gh() are "fine" insofar as they
## contain correct g() and H() formulas

#################################################


##' Cholesky with Hessian =: [ch]o[l]_[H]e[s]sia[n] = chlHsn
##' implementation of the Cholesky decomposition as in nlm()
##'   (see pdf-report Section 2.2)
chlhsn <- function(H, epsm, sx, neworder = FALSE) {
  n <- nrow(H)
  H <- H/(sx %*% t(sx))
  tol <- sqrt(epsm)
  r.d <- range(diag(H)); diagmn <- r.d[1]; diagmx <- r.d[2]
  posmax <- max(0,diagmx)
  if(diagmn <= tol*posmax) { ## FALSE for pos.def. matrices
    amu <- tol*(posmax-diagmn) - diagmn
    if(amu == 0) {
      H1 <- abs(H)
      diag(H1) <- NA
      offmax <- max(0,max(abs(H1,na.rm = TRUE)))
      amu <- if(offmax == 0) 1 else offmax*(1+tol)
    }
    diag(H) <- udiag <- diag(H)+amu
    diagmx <- diagmx+amu
  } else {
    udiag <- diag(H)
    amu <- NA # (to be returned)
  }
  c1 <- choldc(H, diagmx, tol, neworder=neworder) # first attempt of cholesky decomp
  if(c1$addmax > 0) { ## H[] must be modificed to be safely pos.def.
    ## (no need to restore H[] now)
    evmin <- 0
    evmax <- H[1,1]
    for(i in 1:n) { ## Gershgorin circle theorem
      offrow <- sum(abs(H[i,-i]))
      evmin <- min(evmin, H[i,i]-offrow)
      evmax <- max(evmax, H[i,i]+offrow)
    }
    sdd <- tol*(evmax-evmin)-evmin
    amu <- min(sdd, c1$addmax) ## minimal mu to be added to Hessian
    udiag <- (diag(H) <- diag(H) + amu)
    ## cholesky decomposition of modified Hessian:
    c1 <- choldc(H,0,tol, neworder=neworder)
  }
  H <- c1$H
  ## unscale Hessian matrix :
  for(j in 1:n) {
    ii <- j:n
    H[ii,j] <- H[ii,j]*sx[ii]
    if(length(ii <- seq_len(j-1)))
      H[ii,j] <- H[ii,j]*sx[ii]*sx[j]
  }
  udiag <- udiag*sx^2
  list(H = H, udiag = udiag, amu = amu, addmax = c1$addmax)
} ##end {chlhsn}
##
choldc <- function(H, diagmx, tol, neworder = FALSE) {
  addmax <- 0
  aminl <- sqrt(diagmx*tol)
  amnlsq <- aminl^2
  n <- nrow(H)
  for(i in 1:n) {
    jj <- seq_len(i-1) ## j in jj < i
    if(neworder) {
      ## 1) H[i,j] := *,
      for(j in jj) {  ## => j < i
        kk <- seq_len(j-1) # k < j < i
        sumH <- sum(H[i,kk]*H[j,kk])
        H[i,j] <- (H[i,j]-sumH)/H[j,j] # normal Cholesky formula
      }
      ## 2) H[i,i] := *
      tmp1 <- H[i,i] - sum(H[i,jj]^2)
      if(tmp1 >= amnlsq) { # normal Cholesky
        H[i,i] <- sqrt(tmp1)
      } else { ## augment diagonal of H[,]
        offmax <- max(0, abs(H[i,jj]), amnlsq)
        H[i,i] <- sqrt(offmax)
        addmax <- max(addmax, offmax-tmp1)
      }

    } else { ## old _wrong_ order

      sumH <- sum(H[i,jj]^2)
      tmp1 <- H[i,i]-sumH
      if(tmp1 >= amnlsq) { # normal Cholesky
        H[i,i] <- sqrt(tmp1)
      } else { ## augment diagonal of H[,]
        offmax <- max(0, abs(H[i,jj]), amnlsq)
        H[i,i] <- sqrt(offmax)
        addmax <- max(addmax, offmax-tmp1)
      }
      for(j in jj) {  ## => j < i
        kk <- seq_len(j-1) # k < j < i
        sumH <- sum(H[i,kk]*H[j,kk])
        H[i,j] <- (H[i,j]-sumH)/H[j,j] # normal Cholesky formula
      }
    }
  }
  list(H = H, addmax = addmax)
}## {choldc()}
#####################
##
##
## Implementation of the Cholesky decomposition as in Dennis & Schnabel (1996)
##   (see pdf-report Section 2.3)
chlhsnDS <- function(H, epsm, sx, old.chol = FALSE) {
  n <- nrow(H)
  H <- H/(sx %*% t(sx))
  tol <- sqrt(epsm)
  diagmx <- max(diag(H))
  diagmn <- min(diag(H))
  posmax <- max(0,diagmx)
  if(diagmn <= tol*posmax) { ## FALSE for pos.def. matrices
    amu <- 2*tol*(posmax-diagmn) - diagmn
    ##     ==
    diagmx <- diagmx+amu
  } else {
    amu <- 0
  }
  H1 <- abs(H)
  diag(H1) <- NA
  offmax <- max(0, max(H1, na.rm = TRUE))
  if(diagmx <= (offmax*(1+2*tol))) {
    amu <- amu+(offmax-diagmx)+2*tol*offmax
    diagmx <- offmax*(1+2*tol)
  }
  if(diagmx == 0) {
    amu <- 1
    diagmx <- 1
  }
  if(amu > 0) {
    diag(H) <- udiag <- diag(H)+amu
  } else
    udiag <- diag(H)
  maxoffl <- sqrt(max(diagmx, offmax/n))
  c1 <- choldcDS(H, maxoffl, tol)
  if(c1$addmax > 0) {
    ## (no need to restore)
    evmin <- H[1,1]
    evmax <- H[1,1]
    for(i in 1:n) {
      offrow <- sum(abs(H[i,-i]))
      tmp <- H[i,i]-offrow
      evmin <- min(evmin,tmp)
      tmp <- H[i,i]+offrow
      evmax <- max(evmax,tmp)
    }
    sdd <- max(tol*(evmax-evmin)-evmin,0)
    amu <- min(sdd, c1$addmax)
    diag(H) <- diag(H)+amu
    udiag <- diag(H)
    c1 <- choldcDS(H, 0, tol)
  }
  H[lower.tri(H,diag = TRUE)] <- c1$H[lower.tri(H,diag = TRUE)]
  addmax <- c1$addmax
  ## unscale Hessian matrix :
  for(j in 1:n) {
    for(i in j:n) {
      H[i,j] <- H[i,j]*sx[i]
    }
    if(j > 1) {
      for(i in 1:(j-1)) {
        H[i,j] <- H[i,j]*sx[i]*sx[j]
      }
    }
  }
  udiag <- udiag*sx^2
  list(H = H, udiag = udiag, amu = amu, addmax = c1$addmax)
} ## end{ chlhsnDS() }

##' choldc() according to Dennis & Schnabel (1996)
choldcDS <- function(H, maxoffl, tol) {
  L <- array(0, dim=dim(H))
  minl <- maxoffl*sqrt(tol)
  if(maxoffl == 0) maxoffl <- sqrt(max(abs(diag(H))))
  minl2 <- maxoffl*tol
  addmax <- 0
  n <- nrow(H)
  for(j in 1:n) {
    sum <- 0
    if(j > 1) {
      for(i in 1:(j-1)) {
        sum <- sum+L[j,i]^2
      }
    }
    L[j,j] <- H[j,j]-sum
    minljj <- 0
    if(j < n) {
      for(i in (j+1):n) {
        sum2 <- 0
        if(j > 1) sum2 <- sum(L[i,1:(j-1)]*L[j,1:(j-1)])
        L[i,j] <- H[j,i]-sum2
        minljj <- max(minljj, abs(L[i,j]))
      }
    }
    minljj <- max(minljj/maxoffl, minl)
    if(L[j,j] > minljj^2) {
      L[j,j] <- sqrt(L[j,j])
    } else {
      minljj <- max(minljj, minl2)
      addmax <- max(addmax, minljj^2-L[j,j])
      L[j,j] <- minljj
    }
    if(j < n) {
      for(i in (j+1):n) {
        L[i,j] <- L[i,j]/L[j,j]
      }
    }
  }
  list(H = L, addmax = addmax)
}
#####################
##
## application to initial Hessian of Rosenbrock banana valley function
(H <- hessRosenb(-1.2, 1))
eigen(H)$values
kappa(H)
(L <- t(chol(H))) # correct one
epsm <- .Machine$double.eps; tol <- sqrt(epsm)
(diagmx <- max(diag(H)));  diagmn <- min(diag(H));  posmax <- max(0,diagmx)
diagmn <= posmax*tol # FALSE
 choldc(H, diagmx, tol)    ##  wrong [2,2] entry; large addmax = 230680
(choldc(H, diagmx, tol, neworder=TRUE) -> ch2)# correct [2,2] = 5.17; addmax = 0
stopifnot(all.equal(ch2$H[-3], L[-3]), ch2$addmax == 0)
##
## application to initial Hessian of Rosenbrock banana valley function
## with DS [using 'maxoffl' instead of 'diagmx']:
H1 <- abs(H);  diag(H1) <- NA;  offmax <- max(0,max(H1,na.rm = TRUE))
(maxoffl <- sqrt(max(diagmx,offmax/nrow(H)))) # 36.47
(choldcDS(H, maxoffl, tol) -> ch3)
stopifnot(all.equal(ch3$H, L), ch3$addmax == 0)

### MM: application to initial Hessian of Wood function:
(H <- wfh(c(-3,-1,-3,-1)))
eigen(H)$values
kappa(H) # 183.88
(L <- t(chol(H))) # correct one [is easy!]
epsm <- .Machine$double.eps; tol <- sqrt(epsm)
(diagmx <- max(diag(H)));  diagmn <- min(diag(H));  posmax <- max(0,diagmx)
diagmn <= posmax*tol # FALSE
 choldc(H, diagmx, tol)    ##  wrong [2,2] and [4,4] entry; large addmax = 1'440'980
(choldc(H, diagmx, tol, neworder=TRUE) -> ch2) # correct [2,2]=9.573, [4,4]=8.957
lo.tr <- lower.tri(L, diag=TRUE)
stopifnot(all.equal(ch2$H[lo.tr], L[lo.tr])) # correct !!
##
## with DS [using 'maxoffl' instead of 'diagmx']:
H1 <- abs(H);  diag(H1) <- NA;  offmax <- max(0,max(H1,na.rm = TRUE))
(maxoffl <- sqrt(max(diagmx,offmax/nrow(H)))) # 105.8395
(choldcDS(H, maxoffl, tol) -> ch3)
stopifnot(all.equal(ch3$H[lo.tr], L[lo.tr])) # correct, too


##############################################################################
##
## implementation of nlm()-algorithm to check how modified algorithms
##   affect minimization results
##   (see pdf-report Section 3)
##
linesearch <- function(x,f,g,p,fn,state = NULL,stepmx,steptl = 1e-6,sx = 1) {
  firstback <- TRUE
  mxtake <- FALSE
  iretcd <- 2
  sln <- sqrt(sum(sx*p*sx*p))
  if (sln > stepmx) {
    scl <- stepmx/sln
    p <- p*scl
    sln <- stepmx
  }
  slp <- sum(g*p)
  v <- abs(p)/pmax(abs(x),1/sx)
  v <- pmax(v,0)
  rln <- max(v)
  rmnlmb <- steptl/rln
  lambda <- 1
  xpls <- x+lambda*p
  fpls <- fn(xpls)
  ## solution found?
  if (fpls <= f+slp*1e-4*lambda) {
    iretcd <- 0
    if(lambda == 1 && sln > stepmx*0.99) mxtake <- TRUE
    return(list(lambda = lambda, xpls = xpls, fpls = fpls,
                iretcd = iretcd, mxtake = mxtake))
  } else if(lambda < rmnlmb) {
    return(list(lambda = lambda, xpls = xpls, fpls = fpls,
                iretcd = 1, mxtake = mxtake))
  } else {
    if(fpls >= .Machine$double.xmax) {
      lambda <- lambda*0.1
      firstback <- TRUE
    } else {
      ## firstback
      tlmbda <- -lambda*slp/(2*(fpls-f-slp))
      firstback <- FALSE
      plmbda <- lambda
      pfpls <- fpls
      if(tlmbda < (lambda*0.1)) {
        lambda <- lambda*0.1
      } else {
        lambda <- tlmbda
      }
    }
  }
  xpls <- x+lambda*p
  fpls <- fn(xpls)
  while(fpls > f+slp*1e-4*lambda) {
    if (fpls <= f+slp*1e-4*lambda) {
      iretcd <- 0
      if(lambda == 1 && sln > stepmx*0.99) mxtake <- TRUE
      return(list(lambda = lambda,xpls = xpls,fpls = fpls,iretcd = iretcd,
                  mxtake = mxtake))
    }
    if(lambda < rmnlmb)
      return(list(lambda = lambda, xpls = xpls, fpls = fpls,
                  iretcd = 1, mxtake = mxtake))
    if(fpls >= .Machine$double.xmax) {
      lambda <- lambda*0.1
      firstback <- TRUE
    } else {
      if(firstback) {
        tlmbda <- -lambda*slp/(2*(fpls-f-slp))
        firstback <- FALSE
      } else {
        t1 <- fpls-f-lambda*slp
        t2 <- pfpls-f-plmbda*slp
        t3 <- 1/(lambda-plmbda)
        a3 <- 3*t3*(t1/(lambda^2)-t2/(plmbda^2))
        b <- t3*(t2*lambda/(plmbda^2)-t1*plmbda/(lambda^2))
        disc <- b^2-a3*slp
        if(disc > b^2) {
          tlmbda <- (-b+ifelse(a3 < 0,-sqrt(disc),sqrt(disc)))/a3
        } else {
          tlmbda <- (-b+ifelse(a3 < 0,sqrt(disc),-sqrt(disc)))/a3
        }
        if(tlmbda > (0.5*lambda)) tlmbda <- 0.5*lambda
      }
      plmbda <- lambda
      pfpls <- fpls
      if(tlmbda < (0.1*lambda)) {
        lambda <- 0.1*lambda
      } else {
        lambda <- tlmbda
      }
    }
    xpls <- x+lambda*p
    fpls <- fn(xpls)
  }
  list(lambda = lambda,xpls = xpls,fpls = fpls,iretcd = iretcd,mxtake = mxtake)
} ## { linesearch() }
##
opt_stop <- function(xpls,fpls,gpls,x,itncnt,icscmx,gradtl,steptl,sx,fscale,
                     itnlim,iretcd,mxtake,pL) {
  if(iretcd == 1) return(list(itrmcd = 3,icscmx = icscmx))
  d <- max(abs(fpls),fscale)
  rgx <- max(0,max(abs(gpls)*pmax(abs(xpls),1/sx)/d))
  jtrmcd <- 1
  if(rgx > gradtl) {
    if(itncnt == 0) return(list(itrmcd = 0,icscmx = icscmx))
    rsx <- max(0,max(abs(xpls-x)/pmax(abs(xpls),1/sx)))
    jtrmcd <- 2
    if(rsx > steptl) {
      jtrmcd <- 4
      if(itncnt < itnlim) {
        if(!mxtake) {
          return(list(itrmcd = 0, icscmx = 0))
        } else {
          icscmx <- icscmx+1
          if(icscmx < 5) return(list(itrmcd = 0, icscmx = icscmx))
          jtrmcd <- 5
        }
      }
    }
  }
  list(itrmcd = jtrmcd, icscmx = icscmx)
} ## { opt_stop() }
##
optdrv_end <- function(xpls, x, fpls, f, gpls, g, H, p, itncnt, itrmcd, pL) {
  if(itrmcd == 3) {
    fpls <- f
    xpls <- x
    gpls <- g
  }
  if(pL > 0) {
    cat(paste("iteration =",itncnt,"\n"))
    cat("Step:\n");print(p)
    cat("Parameter:\n");print(xpls)
    cat(paste("Function value:\n",fpls,"\n"))
    cat("Gradient:\n");print(gpls)
    if(!is.null(H)) cat("Hessian:\n");print(H)
  }
} ## { optdrv_end() }

##
## without comparison of analytical and numerical gradients and Hessians respectively
optdrv <- function(x, fn, gr = NULL, h = NULL, perturb = FALSE,
                   method = c("DS96", "new.chol", "old.chol"), chk.chlhsn=FALSE,
                   typsiz = rep(1, length(x)), fscale = 1,
                   pL = 0, ndigit = 12, itnlim = 100, gradtl = 1e-6,
                   stepmx = max(1000*sqrt(sum((x/typsiz)^2)), 1000), steptl = 1e-6)
{
  iagflg <- !(is.null(gr))
  ## here:
  if(!is.function(h)) stop("'h' must be a function computing the Hessian matrix")
  iretcd <- 0
  icscmx <- 0
  itncnt <- 0
  fscale <- abs(fscale)
  sx <- abs(1/typsiz)
  mxtake <- FALSE
  epsm <- .Machine$double.eps
  ## skip input checks here...
  n <- length(x)
  p <- rep(0,n)
  rnf <- max(epsm, 10^(-ndigit))#?
  anatl <- max(0.1, sqrt(rnf))
  f <- fn(x)
  if(!iagflg) {
    stop("No gradient provided.")
  } else {
    g <- gr(x)
  }
  iretcd <- -1
  os <- opt_stop(xpls = x, fpls = f, gpls = g, x = x, itncnt, icscmx, gradtl, steptl,
                 sx, fscale, itnlim, iretcd, mxtake, pL)
  itrmcd <- os$itrmcd
  icscmx <- os$icscmx
  if(itrmcd == 0) {
    H <- h(x)
    ## mH <- H
    method <- match.arg(method)
    while(TRUE) {
      itncnt <- itncnt+1
      ## cholesky
      if(perturb) {
        ch <- switch(method,
                     "DS96" = chlhsnDS(H, epsm, sx),
                     "new.chol" = chlhsn(H, epsm, sx, neworder = TRUE),
                     "old.chol" = chlhsn(H, epsm, sx))
        ## H <- ch$H
        L <- ch$H
        if(chk.chlhsn &&
             ## This is triggered indeed, e.g., for
             ##  H <- hessRosenb(-1.022, 1.053); epsm <- 2e-16; sx <- c(1,1)
             !isTRUE(ae <- all.equal(chlhsnDS(H, epsm, sx),
                                     chlhsn  (H, epsm, sx, neworder = TRUE),
                                     tol = 1e-12))) {
          cat("The two version of chlhsn*(H, epsm, sx)  are not equal:\n")
          writeLines(ae)
          browser()# <- you can inspect
        }
        udiag <- ch$udiag
      } else {
        L <- t(chol(H))
      }
      ## solve for Newton step
      y <- forwardsolve(L, -g)
      p <- forwardsolve(L, y, transpose = TRUE)
      ## linesearch
      z <- linesearch(x, f, g, p, fn, state = NULL, stepmx, steptl, sx)
      ## calculate step
      xpls <- z$xpls
      fpls <- z$fpls
      iretcd <- z$iretcd
      mxtake <- z$mxtake
      p <- xpls-x
      gpls <- gr(xpls)
      ## update hessian
      H <- h(xpls)
      ## Ls <- L
      ## Ls[upper.tri(Ls)] <- 0
      ## mH <- Ls %*% t(Ls)
      ## opt_stop
      os <- opt_stop(xpls, fpls, gpls, x = x, itncnt, icscmx, gradtl, steptl,
                     sx, fscale, itnlim, iretcd, mxtake, pL)
      itrmcd <- os$itrmcd
      icscmx <- os$icscmx
      if(itrmcd != 0) break
      ## print
      if(pL > 1) {
        cat("iteration =",itncnt,"\n")
        cat("Step:\n");print(p)
        cat("Parameter:\n");print(xpls)
        cat("Function value:\n",fpls,"\n")
        cat("Gradient:\n");print(gpls)
        cat("Hessian:\n");print(H)
      }
      ## f <- fpls...
      f <- fpls
      x <- xpls
      g <- gpls
    } # end{while}

    optdrv_end(xpls, x, fpls, f, gpls, g, H = H, p, itncnt, itrmcd, pL)
    list(minimum = fpls, estimate = xpls, gradient = gpls, code = itrmcd, iterations = itncnt)

  } else { ## itrmcd != 0:

    optdrv_end(xpls = x, x = NULL, fpls = f, f = NULL, gpls = g, g = NULL, H = NULL,
               p = rep(0, n), itncnt, itrmcd, pL)
    list() # return nothing
  }
} ## end optdrv()

## Rosenbrock banana valley function
f <- function(x) f.Rosenb(x[1], x[2])
## gradient:
gr <- function(x) grRosenb(x[1], x[2])
## hessian:
h <- function(x) hessRosenb(x[1], x[2])
##
## minimization of Rosenbrock banana valley function
res1 <- optdrv(c(-1.2,1), f, gr, h, perturb = TRUE, method="old.chol"); str(res1)
##-> not okay
all.equal(res1, nlm.fgh)# are the same (NOT for R-devel where nlm.fgh is correct!)
res1n<- optdrv(c(-1.2,1), f, gr, h, perturb = TRUE, method="new.chol"); str(res1n)# fine!
res2 <- optdrv(c(-1.2,1), f, gr, h, perturb = TRUE, method = "DS96") ;  str(res2) # fine
all.equal(res1n, res2)# are the same
res3 <- optdrv(c(-1.2,1), f, gr, h, perturb = FALSE)
all.equal(res2,  res3)# are the same  <==>  perturb-chol was not necessary here

#####################
##
## Wood function
##
## minimization of Wood function -- harder: Hessian is _not_ pos.def !
try( reswf0 <- optdrv(c(-3,-1,-3,-1), f.wood, wfg, wfh, perturb=FALSE) )# error in chol()
## ==> here, *do* need 'perturb':
reswf1 <- optdrv(c(-3,-1,-3,-1), f.wood, wfg, wfh, perturb= TRUE, method = "old.chol")
reswf1n<- optdrv(c(-3,-1,-3,-1), f.wood, wfg, wfh, perturb= TRUE, method = "new.chol")
reswf2 <- optdrv(c(-3,-1,-3,-1), f.wood, wfg, wfh, perturb= TRUE, method = "DS96")
str(reswf1) # not ok
str(reswf1n)# _also_ not ok
str(reswf2) #  okay
#####################################################################
