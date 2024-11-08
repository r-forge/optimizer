

require(GNE)
#-------------------------------------------------------------------------------
# (3) River basin pollution game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------

myorigarg <- list(
  C = cbind(c(.1, .12, .15), c(.01, .05, .01)),
  U = cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75)),
  K = c(100, 100),
  E = c(.5, .25, .75),
  D = c(3, .01)
)

#### objective functions ####

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

#### contraint functions ####
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

#### GNE computation ####

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

if(FALSE)
{
  
  basesetting <- matrix(unlist(myorigarg), nrow=1, dimnames = list(NULL, names(unlist(myorigarg))))
  resbase <- getNE(basesetting, output="full")
  
  resprofit <- sapply(1:3, function(j) 
    obj(resbase[,1]$par, j, myorigarg))
  
  resconstr <- h(resbase[,1]$par, myorigarg)
  
}

#### reverse engineering for nu-beta ####

#sj= position of agent j
#sl= position of quality control l
expdist <- function(nu, sj, sL)
{
  #Equation (11)
  #A. Haurie, J.B. Krawczyk / Optimal charges on river effluent 
  #set beta=1
  exp(-nu*(sL-sj))*(sL >= sj)
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
    T11 <- alphatab[1,1] - expdist(nu,  sJ1, sL1)
    T12 <- alphatab[1,2] - expdist(nu,  sJ1, sL2)
    T21 <- alphatab[2,1] - expdist(nu,  sJ2, sL1)
    T22 <- alphatab[2,2] - expdist(nu,  sJ2, sL2)
    T31 <- alphatab[3,1] - expdist(nu,  sJ3, sL1)
    T32 <- alphatab[3,2] - expdist(nu,  sJ3, sL2)
    
    (T11^2+T12^2+T21^2+T22^2+T31^2+T32^2)/2
  }
  gradf <- function(x)
  {
    nu <- x[1]
    beta <- x[2]
    sJ1 <- x[3]
    sJ2 <- x[4]
    sJ3 <- x[5]
    sL1 <- x[6]
    sL2 <- x[7]
    alpha11 <- expdist(nu,  sJ1, sL1)
    alpha12 <- expdist(nu,  sJ1, sL2)
    alpha21 <- expdist(nu,  sJ2, sL1)
    alpha22 <- expdist(nu,  sJ2, sL2)
    alpha31 <- expdist(nu,  sJ3, sL1)
    alpha32 <- expdist(nu,  sJ3, sL2)
    T11 <- alphatab[1,1] - alpha11
    T12 <- alphatab[1,2] - alpha12
    T21 <- alphatab[2,1] - alpha21
    T22 <- alphatab[2,2] - alpha22
    T31 <- alphatab[3,1] - alpha31
    T32 <- alphatab[3,2] - alpha32
    
    M <- rbind(
      c(sL1-sJ1, sJ1-sL2, sJ2-sL1, sJ2-sL2, sJ3-sL1, sJ3-sL2)/beta, #w.r.t. nu
      -c(sL1-sJ1, sJ1-sL2, sJ2-sL1, sJ2-sL2, sJ3-sL1, sJ3-sL2)/beta^2, #w.r.t. beta
      c(-nu/beta, -nu/beta, 0, 0, 0, 0), #w.r.t. sJ1
      c(0, 0, -nu/beta, -nu/beta, 0, 0), #w.r.t. sJ2
      c(0, 0, 0, 0, -nu/beta, -nu/beta), #w.r.t. sJ3
      c(nu/beta, 0, nu/beta, 0, nu/beta, 0), #w.r.t. sL1
      c(0, nu/beta, 0, nu/beta, 0, nu/beta) #w.r.t. sL2
    )
    M <- M %*% diag(c(T11, T12, T21, T22, T31, T32))
    M <- M %*% diag(c(alpha11, alpha12, alpha21, alpha22, alpha31, alpha32))
    M %*% rep(1, 6) #sum the six terms
  }
  x0 <- c("nu"=1, "beta"=1, "sJ"=1:3/4, "sL"=7:8/8)
  cat(f(x0), "\n")
  print(gradf(x0))
  
  matU <- rbind(
    c(1, rep(0, 6)),
    c(0, 1, rep(0, 5)),
    cbind(0, 0, rep(-1, 2), rep(0, 2), rep(0, 2), diag(2)),
    cbind(0, 0, rep(0, 2), rep(-1, 2), rep(0, 2), diag(2)),
    cbind(0, 0, rep(0, 2), rep(0, 2), rep(-1, 2), diag(2))
  )
  
  constrOptim(x0, f, grad=gradf, ui=matU, ci=rep(0, NROW(matU)), ...)
}

#nubetaNM <- getnubeta(myorigarg$U, method="Nelder-Mead", control=list(trace=1, REPORT=1))
nubetaBFGS <- getnubeta(myorigarg$U, method="BFGS", control=list(trace=1, REPORT=1))

nubetaBFGS$par


#### parameter random generators ####


# parameter U_j,l
gen_Ujl <- function(n, river=list(begin=0, end=1, quality=2/3), nbagent, nbcontrol,
                    nu, beta, echo=FALSE)
{
  stopifnot(is.list(river))
  stopifnot(all(names(river) %in% c("begin", "end", "quality")))
  stopifnot(river$quality < river$end && 0 < river$quality)
  stopifnot(is.numeric(nbagent))
  stopifnot(is.numeric(nbcontrol))
  
  agentpos <- matrix(runif(nbagent * n, min=river$begin, max=river$quality), 
                     nrow=n, ncol=nbagent)
  qualitypos <- matrix(runif(nbcontrol * n, min=river$quality, max=river$end), 
                       nrow=n, ncol=nbcontrol)
  
  if(echo)
  {
    print(cbind.data.frame("agent"=agentpos, "control"=qualitypos))
  }
  j <- l <- 1
  res <- matrix(NA, nrow=n, ncol=nbagent*nbcontrol,
                dimnames=list(1:n, paste0("U", 1:nbagent, rep(1:nbcontrol, each=nbagent))))  
  for(i in 1:n)
  {
    for(l in 1:nbcontrol)
    {
      ulj <- sapply(1:nbagent, function(j) expdist(nu,  agentpos[i,j], qualitypos[i,l]) )
      res[i,1:nbagent + (l-1)*nbagent] <- ulj
      
    }
  }
  res
}

gen_Ujl(5, nbagent=3, nbcontrol=2, nu = 1.05,
        echo=FALSE)





par(mfrow=c(1,1), mar=c(4,4,2,1))

gend1d2 <- function(n, d1.low, d1.upp, d2.low, d2.upp)
{
  res <- cbind(runif(n, d1.low, d1.upp), runif(n, d2.low, d2.upp))
  colnames(res) <- c("D1", "D2")
  res
}
pairs(gend1d2(1000, 2, 4, 0.01, 0.02))

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
gencij(3, as.vector(myorigarg$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01))
pairs(gencij(1e3, as.vector(myorigarg$C), 0, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01)))

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
getUj(10, as.vector(myorigarg$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01))

pairs(getUj(1e3, as.vector(myorigarg$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01)))

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

pairs(getKj(1e3, myorigarg$K, 0.5, c(5, 5)))

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

cstU <- as.vector(myorigarg$U)
names(cstU) <- paste0("U", 1:6)

#test
X1 <- genXparam(10, as.vector(myorigarg$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                as.vector(myorigarg$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                myorigarg$K, 0.5, c(5, 5), .1, .9,
                myorigarg$D, 0.5, c(0.10, 0.025)) 

fX1 <- getNE(X1[1:10, ])
dim(fX1)

n <- 1e2
X1 <- genXparam(n, as.vector(myorigarg$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                as.vector(myorigarg$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                myorigarg$K, 0.5, c(5, 5), .1, .9,
                myorigarg$D, 0.5, c(0.10, 0.025)) 
dim(X1)
X2 <- genXparam(n, as.vector(myorigarg$C), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                as.vector(myorigarg$U), 0.5, c(0.1, 0.1, 0.1, 0.01, 0.01, 0.01),
                myorigarg$K, 0.5, c(5, 5), .1, .9,
                myorigarg$D, 0.5, c(0.10, 0.025)) 
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
    obj(Youtputs[i,], j, myorigarg)
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

