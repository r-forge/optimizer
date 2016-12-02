## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----simpleproduct, echo=TRUE, cache=TRUE--------------------------------
cat("try loading optimrx\n")
require(optimrx, quietly=TRUE)
pr <- function(y) {
- prod(y)*(1-sum(y))
}
cat("test the simple product for n=5\n")
meth <- c("Nelder-Mead", "BFGS")
n<-5
  st<-1:(n-1)/(n*n)
   ans<-opm(st, pr, gr="grcentral", control=list(trace=0))
   ao<-summary(ans,order=value)
print(ao)

## ----gabornll, echo=TRUE, cache=TRUE-------------------------------------
nll <- function(y) {
  if ((any(y <= 10*.Machine$double.xmin)) || (sum(y)>1-.Machine$double.eps))
         .Machine$double.xmax
  else   - sum(log(y)) - log(1-sum(y))
}
nll.g <- function(y) { - 1/y + 1/(1-sum(y))} # so far not safeguarded

## ----label=C13badruns1, echo=TRUE, cache=TRUE----------------------------
require(optimrx, quietly=TRUE)
n<-5
mset<-c("L-BFGS-B", "BFGS", "CG", "spg", "ucminf", "nlm", "nlminb", "Rvmmin", "Rcgmin")
a5<-opm(2:n/n^2, nll, gr="grfwd", method=mset, control=list(dowarn=FALSE))
a5g<-opm(2:n/n^2, nll, nll.g, method=mset, control=list(dowarn=FALSE))
a5gb<-opm(2:n/n^2, nll, nll.g, lower=0, upper=1, method=mset, control=list(dowarn=FALSE))
#- a5x <- opm(2:n/n^2, nll, nll.g, method="ALL", control=list(dowarn=FALSE))
summary(a5,order=value)
summary(a5g,order=value)
summary(a5gb,order=value)
#- summary(a5x,order=value)

## ----label=C13ravi1, echo=TRUE, cache=TRUE-------------------------------
require(BB, quietly=TRUE)
nllrv <- function(x) {- sum(log(x))}
nllrv.g <- function(x) {- 1/x }
proj <- function(x) {x/sum(x)}
n <- 5
tspg<-system.time(aspg <- spg(par=(1:n)/n^2, fn=nllrv, gr=nllrv.g, project=proj))[[3]]
tspgn<-system.time(aspgn <- spg(par=(1:n)/n^2, fn=nllrv, project=proj))[[3]]
cat("Times: with gradient =",tspg,"   using numerical approx.=", tspgn,"\n")
cat("F_optimal: with gradient=",aspg$value,"  num. approx.=",aspgn$value,"\n")
pbest<-rep(1/n, n)
cat("fbest = ",nllrv(pbest),"  when all parameters = ", pbest[1],"\n")
cat("deviations:  with gradient=",max(abs(aspg$par-pbest)),"   num. approx.=",max(abs(aspg$par-pbest)),"\n")

## ----expgabor, echo=TRUE, cache=TRUE-------------------------------------
enll <- function(lx) {
    x<-exp(lx)
    fval<-  - sum( log( x/sum(x) ) ) 
}
enll.g <- function(lx){
    x<-exp(lx)
    g<-length(x)/sum(x) - 1/x
    gval<-g*exp(lx)
}

## ----expgabrun1, warning=FALSE, echo=TRUE, cache=TRUE--------------------
require(optimrx, quietly=TRUE) # just to be sure
st<-1:5/10 # 5 parameters, crude scaling to start
a5x<-opm(st, enll, enll.g, method="ALL", control=list(trace=0))
a5xbyvalue<-summary(a5x, order=value)
xnor<-a5xbyvalue[1, 1:5] # get the 5 parameters of "best" solution
xnor<-xnor/sum(xnor)
cat("normalized parameters:")
print(xnor)

## ----expgabrun2, warning=FALSE, echo=TRUE, cache=TRUE--------------------
require(Rcgmin, quietly=TRUE)
st<-1:100/1e3 # large
stenll<-enll(st)
cat("Initial function value =",stenll,"\n")
tym<-system.time(acgbig<-Rcgmin(st, enll, enll.g, control=list(trace=0, tol=1e-32)))[[3]]
cat("Time = ",tym,"  fval=",acgbig$value,"\n")
xnor<-acgbig$par/sum(acgbig$par)
print(xnor)

## ----sphere0, echo=FALSE, cache=TRUE-------------------------------------
library(BB)
library(optimx)

## ----sphere5, echo=TRUE, cache=TRUE--------------------------------------
proj2 <- function(theta) {
    theta2 <- theta^2
    s2 <- theta2 / (1 + theta2)
    cumprod(c(1, s2)) * c(1-s2, 1)
 }
obj <- function(theta) - sum(log(proj2(theta)))
 n <- 5
 ans <- spg(seq(n-1), obj)
 proj2(ans$par)

## ----sphere100, echo=TRUE, cache=TRUE------------------------------------
n<-100
ans100 <- spg(seq(n-1), obj, control=list(trace=FALSE), quiet=TRUE)
proj2( (ans100$par) )

## ----sphere100all, echo=TRUE, cache=TRUE---------------------------------
allans<- opm(seq(n-1), obj, gr="grfwd", method="ALL", control=list(dowarn=FALSE))
summary(allans, order = "list(round(value, 3), fevals)", par.select = FALSE)

## ----rayspg1, echo=TRUE, cache=TRUE--------------------------------------
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

raynum<-function(x, A){
   rayquo<-as.numeric((t(x)%*%A)%*%x)
}

proj<-function(x) { x/sqrt(crossprod(x)) }

require(BB, quietly=TRUE)
n<-10
x<-rep(1,n)
A<-molerbuild(n)
tmin<-system.time(asprqmin<-spg(x, fn=raynum, project=proj, A=A))[[3]]
tmax<-system.time(asprqmax<-spg(x, fn=raynum, project=proj, A=-A))[[3]]
cat("maximal eigensolution: Value=",asprqmax$value,"in time ",tmax,"\n")
print(asprqmax$par)
cat("minimal eigensolution: Value=",asprqmin$value,"in time ",tmin,"\n")
print(asprqmin$par)

## ----tran1, echo=TRUE, cache=TRUE----------------------------------------
ssums<-function(x){
  n<-length(x)
  tt<-sum(x)
  ss<-1:n
  xx<-(x/tt)*(x/tt)
  sum(ss*xx)
}

cat("Try penalized sum\n")
require(optimx)
st<-runif(3)
aos<-opm(st, ssums, gr="grcentral", method="ALL")
# rescale the parameters
nsol<-dim(aos)[1]
for (i in 1:nsol){ 
  tpar<-aos[i,1:3] 
  ntpar<-sum(tpar)
  tpar<-tpar/ntpar
#  cat("Method ",aos[i, "meth"]," gives fval =", ssums(tpar))
  aos[i, 1:3]<-tpar 
  
}

summary(aos,order=value)[1:5,]

## ----transpg1, echo=TRUE, cache=TRUE-------------------------------------
ssum<-function(x){
  n<-length(x)
  ss<-1:n
  xx<-x*x
  sum(ss*xx)
}
proj.simplex <- function(y) {
# project an n-dim vector y to the simplex Dn
# Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
# Ravi Varadhan, Johns Hopkins University
# August 8, 2012

n <- length(y)
sy <- sort(y, decreasing=TRUE)
csy <- cumsum(sy)
rho <- max(which(sy > (csy - 1)/(1:n)))
theta <- (csy[rho] - 1) / rho
return(pmax(0, y - theta))
}
as<-spg(st, ssum, project=proj.simplex)
cat("Using project.simplex with spg: fmin=",as$value," at \n")
print(as$par)

## ----label=TranQP, echo=TRUE, cache=TRUE---------------------------------
library(quadprog)
Dmat<-diag(c(1,2,3))
Amat<-matrix(c(1, 1, 1), ncol=1)
bvec<-c(1)
meq=1
dvec<-c(0, 0, 0)
ans<-solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
ans

