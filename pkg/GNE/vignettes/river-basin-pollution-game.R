

require(GNE)
#-------------------------------------------------------------------------------
# (3) River basin pollution game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------

myarg0 <- list(
  C = cbind(c(.1, .12, .15), c(.01, .05, .01)),
  U = cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75)),
  K = c(100, 100),
  E = c(.5, .25, .75),
  D = c(3, .01)
)

dimx <- c(1, 1, 1)
#O_i(x)
obj <- function(x, j, arg)
{
  (arg$D[1] - arg$D[2]*sum(x[1:3]) - arg$C[j, 1] - arg$C[j, 2]*x[j]) * x[j]
}
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
  dij <- 1*(i == j)
  res <- -(-arg$D[2] - arg$C[i, 2]*dij) * x[i] 
  res - (arg$D[1] - arg$D[2]*sum(x[1:3]) - arg$C[i, 1] - arg$C[i, 2]*x[i]) * dij
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
  dij <- 1*(i == j)
  dik <- 1*(i == k)
  
  arg$D[2] * dik + arg$D[2] * dij + 2 * arg$C[i, 2] * dij * dik
}

dimmu <- 5
  #h(x)
  h <- function(x, arg)
    c(sum(arg$U[, 1] * arg$E * x[1:3]) - arg$K[1],
      sum(arg$U[, 2] * arg$E * x[1:3]) - arg$K[2],
      -x[1],
      -x[2],
      -x[3])
  #Gr_x_j h(x)
  grh <- function(x, j, arg)
    c(arg$U[j, 1] * arg$E[j], 
      arg$U[j, 2] * arg$E[j], 
      -1*(1 ==j), 
      -1*(2 ==j), 
      -1*(3 ==j))
  #Gr_x_k Gr_x_j g_i(x)
  heh <- function(x, j, k, arg)
    c(0, 0, 0, 0, 0)
  

#true value around (21.146, 16.027, 2.724, 0.574, 0.000)
z0 <- rep(1, sum(dimx)+sum(dimmu))

getNE <- function(x, control=list(maxit=100, trace=0), check=TRUE, output="equilibrium")
{
  res <- sapply(1:NROW(x), function(i)
    {
      myarg <- list(
        C = cbind(x[i,paste("C",1:3,sep="")], x[i,paste("C",1:3+3,sep="")]),
        U = cbind(x[i,paste("U",1:3,sep="")], x[i,paste("U",1:3+3,sep="")]),
        K = x[i,paste("K",1:2,sep="")],
        E = x[i,paste("E",1:3,sep="")],
        D = x[i,paste("D",1:2,sep="")]
      )
      z0 <- rep(1, sum(dimx)+sum(dimmu))
      
      res <- GNE.nseq(z0, dimx, dimmu=dimmu, grobj=grobj, arggrobj=myarg, heobj=heobj, argheobj=myarg, 
                      joint=h, argjoint=myarg, grjoint=grh, arggrjoint=myarg, hejoint=heh, arghejoint=myarg, 
                      compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
                      control=control)
      
      if(output == "equilibrium")
        return(res$par[1:3])
      else
        return(res)
    }
    )
    t(res)
}

basesetting <- matrix(unlist(myarg0), nrow=1, dimnames = list(NULL, names(unlist(myarg0))))
resbase <- getNE(basesetting, output="full")

resprofit <- sapply(1:3, function(j) 
  obj(resbase[,1]$par, j, myarg0))

resconstr <- h(resbase[,1]$par, myarg0)


myarg0


expdist <- function(nu, beta, sj, sL)
{
  #equation (11)
  exp(-nu/beta*(sL-sj))*(sL >= sj)
}

getnubeta <- function(alphatab, ...)
{
  
  f <- function(x)
  {
    nu <- x[1]
    beta <- x[2]
    sJ1 <- x[3]
    sJ2 <- x[4]
    sJ3 <- x[5]
    sL1 <- x[6]
    sL2 <- x[7]
    T11 <- alphatab[1,1] - expdist(nu, beta, sJ1, sL1)
    T12 <- alphatab[1,2] - expdist(nu, beta, sJ1, sL2)
    T21 <- alphatab[2,1] - expdist(nu, beta, sJ2, sL1)
    T22 <- alphatab[2,2] - expdist(nu, beta, sJ2, sL2)
    T31 <- alphatab[3,1] - expdist(nu, beta, sJ3, sL1)
    T32 <- alphatab[3,2] - expdist(nu, beta, sJ3, sL2)
    
    T11^2+T12^2+T21^2+T22^2+T31^2+T32^2
  }
  x0 <- c("nu"=10, "beta"=10, "sJ"=1:3, "sL"=7:8)
  cat(f(x0), "\n")
  optim(x0, f, ...)
}

nubetaNM <- getnubeta(myarg0$U)
nubetaBFGS <- getnubeta(myarg0$U, method="BFGS")

nubetaBFGS <- getnubeta(myarg0$U, lower=c(1, 1, rep(0,3), 1, 1), 
                        method="L-BFGS-B", control=list(trace=1, REPORT=1))


cbind(nubetaNM$par, nubetaBFGS$par)



dStationJoueur1 <- function(x) 
  expdist(nubetaBFGS$par["nu"], nubetaBFGS$par["beta"], nubetaBFGS$par["sJ1"], x)
dStationJoueur2 <- function(x) 
  expdist(nubetaBFGS$par["nu"], nubetaBFGS$par["beta"], nubetaBFGS$par["sJ2"], x)

curve(dStationJoueur1(x), from = 0, to = 10, n =501, xlab=expression(s),
      ylab=expression(alpha[list(j,l)](s)))
curve(dStationJoueur2(x), add=TRUE, col="blue", n=501)
abline(v=nubetaBFGS$par[3:4], lty=2)


qbar <- function(sJ1, sJ2, sJ3, sL, emrate, prodlevel, nu, beta)
{
  alpha1L <- expdist(nu, beta, sJ1, sL)
  alpha2L <- expdist(nu, beta, sJ2, sL)
  alpha3L <- expdist(nu, beta, sJ3, sL)
  #cat(emrate, "\n")
  #cat(prodlevel, "\n")
  resterms <- c(alpha1L, alpha2L, alpha3L) * emrate * prodlevel
  #cat(resterms, "\n")
  sum(resterms)
}


qbarf <- function(x, resopt, emrate, prodlevel)
{
  sapply(x, function(pos)
    qbar(resopt["sJ1"], resopt["sJ2"], resopt["sJ3"], pos, 
       emrate, prodlevel, resopt["nu"], resopt["beta"])
  )
}

qbarf(3, nubetaBFGS$par, basesetting[1, c("E1", "E2", "E3")], resbase[1,]$par[1:3])

#resbase[1,]$par[1:3]
#basesetting[1, c("E1", "E2", "E3")],

curve(qbarf(x, nubetaBFGS$par, c(1,1,1),  c(1,1,1)),
      xlab=expression(s), ylab=expression(bar(q)(s)), from=0, to=10, n=1001)
abline(v=nubetaBFGS$par[3:5], lty=2, col="grey")
abline(h=100, col="red")

curve(qbarf(x, nubetaNM$par, c(.5, .25, .75), resbase[1,]$par[1:3]),
      from=0, to=10, n=1001)
abline(v=nubetaNM$par[3:5], lty=2, col="grey")
abline(h=100, col="red")







par(mfrow=c(1,1), mar=c(4,4,2,1))

gend1d2 <- function(n, d1d2.mean, d1d2.rho, d1d2.sig)
{
  require(mvtnorm)
  Sig <- diag(d1d2.sig) %*% cbind(c(1, d1d2.rho), c(d1d2.rho, 1)) %*% diag(d1d2.sig)
  #print(Sig)
  res <- rmvnorm(n, d1d2.mean, Sig)
  #print(res)
  res <- pmax(res, 0.005)
  colnames(res) <- c("D1", "D2")
  res
}
pairs(gend1d2(1000, myarg0$D, 0.5, c(0.10, 0.025)))

gencij <- function(n, cij.mean, cij.rho, cij.sig)
{
  require(mvtnorm)
  R <- matrix(cij.rho, 6, 6)
  diag(R) <- 1
  #print(R)
  Sig <- diag(cij.sig) %*% R %*% diag(cij.sig)
  #print(Sig)
  res <- rmvnorm(n, cij.mean, Sig)
  #print(res)
  res <- pmax(res, 0.005)
  #lapply(1:NROW(res), function(i)
  #{
  #  xi <- matrix(res[i,], nrow=3)
  #  colnames(xi) <- c("c1j", "c2j")
  #  xi
  #})
  colnames(res) <- paste0("C", 1:6)
  res
}
gencij(3, as.vector(myarg0$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01))
pairs(gencij(1e3, as.vector(myarg0$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01)))

getej <- function(n, low, upp)
{
  res <- cbind(runif(n, low, upp), runif(n, low, upp), runif(n, low, upp))
  colnames(res) <- paste0("E", 1:3)
  res
}
pairs(getej(1e3, .1, .9))


getUj <- function(n, Uj.mean, Uj.rho, Uj.sig)
{
  require(mvtnorm)
  R <- matrix(Uj.rho, 6, 6)
  diag(R) <- 1
  #print(R)
  Sig <- diag(Uj.sig) %*% R %*% diag(Uj.sig)
  #print(Sig)
  res <- rmvnorm(n, Uj.mean, Sig)
  #print(res)
  res <- pmax(res, 0.005)
  colnames(res) <- paste0("U", 1:6)
  res
}

pairs(getUj(1e3, as.vector(myarg0$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01)))

getKj <- function(n, Kj.mean, Kj.rho, Kj.sig)
{
  require(mvtnorm)
  Sig <- diag(Kj.sig) %*% cbind(c(1, Kj.rho), c(Kj.rho, 1)) %*% diag(Kj.sig)
  #print(Sig)
  res <- rmvnorm(n, Kj.mean, Sig)
  #print(res)
  res <- pmax(res, 0.005)
  colnames(res) <- c("K1", "K2")
  res
}

pairs(getKj(1e3, myarg0$K, 0.5, c(5, 5)))

genXparam <- function(n, cij.mean, cij.rho, cij.sig, Uj.mean, Uj.rho, Uj.sig, 
                      Kj.mean, Kj.rho, Kj.sig, 
                      low, upp, d1d2.mean, d1d2.rho, d1d2.sig)
{
  P1 <- gencij(n, cij.mean, cij.rho, cij.sig)
  P2 <- getUj(n, Uj.mean, Uj.rho, Uj.sig)
  P3 <- getKj(n, Kj.mean, Kj.rho, Kj.sig)
  P4 <- getej(n, low, upp)
  P5 <- gend1d2(n, d1d2.mean, d1d2.rho, d1d2.sig)
  cbind(P1, P2, P3, P4, P5)
}

cstU <- as.vector(myarg0$U)
names(cstU) <- paste0("U", 1:6)

#test
X1 <- genXparam(10, as.vector(myarg0$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          as.vector(myarg0$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          myarg0$K, 0.5, c(5, 5), .1, .9,
          myarg0$D, 0.5, c(0.10, 0.025)) 

fX1 <- getNE(X1[1:10, ])
dim(fX1)

n <- 1e2
X1 <- genXparam(n, as.vector(myarg0$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          as.vector(myarg0$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          myarg0$K, 0.5, c(5, 5), .1, .9,
          myarg0$D, 0.5, c(0.10, 0.025)) 
dim(X1)
X2 <- genXparam(n, as.vector(myarg0$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          as.vector(myarg0$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
          myarg0$K, 0.5, c(5, 5), .1, .9,
          myarg0$D, 0.5, c(0.10, 0.025)) 
dim(X2)


library(sensitivity)
multSob1 <- sobolMultOut(getNE, q=3, X1 = X1, X2 = X2,
                         MCmethod = "sobol")
multSob1
plot(multSob1, ylim=range(multSob1$S))
grid()


Youtputs <- getNE(X1)

ObjOutputs <- sapply(1:3, function(j)
  sapply(1:NROW(Youtputs), function(i)
  obj(Youtputs[i,], j, myarg0)
  )
)

Objarithmean <- apply(ObjOutputs, 1, mean)

(respccArith <- pcc(data.frame(X1), Objarithmean))

(ressrcArith <- src(data.frame(X1), Objarithmean))

plot(respccArith)
grid()
plot(ressrcArith)
grid()




### Arithmetic mean of equilibrium $\bx^\star$

Yarithmean <- apply(Youtputs, 1, mean)

(respccArith <- pcc(data.frame(X1), Yarithmean))

(ressrcArith <- src(data.frame(X1), Yarithmean))

plot(respccArith)
grid()
plot(ressrcArith)
grid()


### Geometric mean of equilibrium $\bx^\star$


Ygeomean <- apply(Youtputs, 1, function(x) prod(pmax(x,0)^(1/3)))

(respccGeo <- pcc(data.frame(X1), Ygeomean))

(ressrcGeo <- src(data.frame(X1), Ygeomean))
  
plot(respccGeo)
grid()

