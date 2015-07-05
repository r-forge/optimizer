projectLinear <- function(par, A, b, meq) {
# A projection function to incorporate linear equalities and 
# inequalities in nonlinear optimization using `spg'
#
# The inequalities are defined such that:  A %*% x - b > 0

   # local function
   solve.QP <- function(dvec, Amat, bvec, meq=0, EPS=1e-07){
     sol <- dvec
     imeq <- seq_len(meq)

     Nmat <- NULL
     wvec <- NULL
     active <- NULL

     repeat{
       viol <-  crossprod(Amat, sol) - bvec
       viol1 <- viol / pmax(1,abs(bvec))
       ## viol1 <- viol / (1 + abs(bvec))

       iim <- viol1[imeq] >= EPS
       if( any(iim) ){
   	iim <- which(iim)
   	viol[iim] <- -viol[iim]
   	bvec[iim] <- -bvec[iim]
   	Amat[,iim] <- -Amat[,iim]
       }

       ii <- which.min(viol1)[1]
       if( viol1[ii] > -EPS) break

       if(ii %in% active)
   	stop("Error in projection")

       wvec <- c(wvec, 0)
       active <- c(active, ii)
       npvec <- Amat[,ii]
       if( !is.null(Nmat) ){
   	rvec <- solve(qr(Nmat, LAPACK=TRUE), npvec)
   	dvec <- npvec - c(Nmat %*% rvec)
       }else{
   	dvec <- npvec
   	rvec <- NULL
       }

       jj <- rvec > 0
       jj[1] <- FALSE

       tmp <- wvec[jj]/rvec[jj]
       t1 <- suppressWarnings(min(tmp))
       t2 <- -viol[ii]/crossprod(npvec, dvec)

       t <- min(c(t1, t2))
       if( !is.finite(t) || t < 0 || t1 <= t2 )
   	stop("Error in projection")

       sol <- sol + t * dvec
       wvec <- wvec - t * c(rvec, -1)
       Nmat <- cbind(Nmat, npvec)
     }
     return(c(sol))
   }

n <- length(par)
if (meq > 0 | any(b - c(A %*% par) > 0)){
    ans <-solve.QP(dvec=rep(0, n), Amat=t(A),
                       bvec = b - c(A %*% par), meq=meq, EPS=1e-07)
    par <- par + ans
}
par
}
