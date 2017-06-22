source("gersub.R", echo=TRUE)
library(geigen)
## rungeradin1
molerbuild<-function(n){ # Create the moler matrix of order n
  # A[i,j] = i for i=j, min(i,j)-2 otherwise
  A <- matrix(0, nrow = n, ncol = n)
  j <- 1:n
  for (i in 1:n) {
    A[i, 1:i] <- pmin(i, 1:i) - 2
  }
  A <- A + t(A)
  diag(A) <- 1:n
  A
}

frankbuild<-function(n){ # Create the Frank matrix of order n
  # A[i,j] = min(i,j)
  cat("Build frank matrix of order",n,"\n")
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
     for (j in 1:n){ A[i,j] <- min(i, j)}
  }
  A
}

unitbuild <- function(n){ A <- diag(rep(1,n)) }

geigres <- function(AA, BB, x, eval){
    # assume x standardized
    res <- (AA - eval * BB) %*% x
    max(abs(res))
}

stdvec <- function(x) {
  # standardize an eigenvector
  vsign <- 1
  bigx <- max(abs(x))
  ix <- which(abs(x) == bigx)
  if (x[ix] < 0) vsign <- -1
  x <- vsign * x /sqrt(as.numeric(crossprod(x)))
  x  
}

cmpesol <- function(x, xcmp, e, ecmp){
  # compare eigensolution to a standard one
  x <- stdvec(x)
  xcmp <- stdvec(xcmp)
  ediff <- abs(e-ecmp)
  vdiff <- max(abs(x - xcmp))
  cat("Comparison: eigenvalue diff=",ediff,"  vector max diff=",vdiff,"\n")
  list(ediff, vdiff)
}

cat("Test geradin eigenvalue method\n")
n <- as.numeric(readline("order=?"))
# ?? test integer, positive?

cat("Frank matrix\n")
AA <- frankbuild(n)
BB=unitbuild(n)

x<-stdvec(rep(1,n) + 0.1*runif(n))

#cat("B Matrix:\n")
#print(BB)
cat("xstart:")
print(x)
tg<-system.time(ag<-geradin(x, ax, bx, AA=AA, BB=BB, 
   control=list(trace=TRUE)))[[3]]
cat("Minimal eigensolution\n")
ag$x<-stdvec(ag$x) # rescale
print(ag)
cat("Geradin time=",tg,"Eigenvalue=",ag$RQ,"\n")
te<-system.time(esol<-geigen(AA,BB, symmetric=TRUE))
cat("geigen time for all eigensolutions:",te[[3]],"\n")
esol <- eigen(AA)
esolmaxe <- esol$values[1]
esolmaxv <- stdvec(esol$vectors[,1])
esolmine <- esol$values[n]
esolminv <- stdvec(esol$vectors[,n])
cat("Maxres - Geradin:", geigres(AA, BB, ag$x, ag$RQ),"   geigen:",
           geigres(AA, BB, esolminv, esolmine),"\n")

tgn<-system.time(agn<-geradin(x, ax, bx, AA=-AA, BB=BB,
   control=list(trace=TRUE)))[[3]]
cat("Maximal eigensolution (negative matrix)\n")
agn$x<-stdvec(agn$x) # rescaleA
print(agn)
cat("Geradin time=",tgn,"Eigenvalue=",-agn$RQ,"\n")
cat("Maxres - Geradin:", geigres(AA, BB, agn$x, -agn$RQ),"   geigen:",
    geigres(AA, BB, esolmaxv, esolmaxe),"\n")

-- more or less ok to here --
cat("Maximal solution - value=",esol$values[1],"  diff =",esol$values[1]+agn$RQ,"\n")
# Note using negative matrix, so use +
print(esol$vectors[,1])
cat("max(abs(diff))=", max(abs(esol$vectors[,1]-agn$x)),"\n")
cat("Minimal solution - value=",esol$values[n],"  diff =",esol$values[n]-ag$RQ,"\n")
# Note using negative matrix, so use +
print(esol$vectors[,n])
cat("max(abs(diff))=", max(abs(esol$vectors[,n]-ag$x)),"\n")
