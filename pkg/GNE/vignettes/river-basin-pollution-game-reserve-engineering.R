

myorigarg <- list(
  C = cbind(c(.1, .12, .15), c(.01, .05, .01)),
  U = cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75)),
  K = c(100, 100),
  E = c(.5, .25, .75),
  D = c(3, .01)
)

#### reverse engineering for nu-beta ####

expdist <- function(nu, beta, sj, sL)
{
  #equation (11)
  exp(-nu/beta*(sL-sj))*(sL >= sj)
}
expdist2 <- function(nu, sj, sL)
{
  #equation (11)
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
    T11 <- alphatab[1,1] - expdist(nu, beta, sJ1, sL1)
    T12 <- alphatab[1,2] - expdist(nu, beta, sJ1, sL2)
    T21 <- alphatab[2,1] - expdist(nu, beta, sJ2, sL1)
    T22 <- alphatab[2,2] - expdist(nu, beta, sJ2, sL2)
    T31 <- alphatab[3,1] - expdist(nu, beta, sJ3, sL1)
    T32 <- alphatab[3,2] - expdist(nu, beta, sJ3, sL2)
    
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
    alpha11 <- expdist(nu, beta, sJ1, sL1)
    alpha12 <- expdist(nu, beta, sJ1, sL2)
    alpha21 <- expdist(nu, beta, sJ2, sL1)
    alpha22 <- expdist(nu, beta, sJ2, sL2)
    alpha31 <- expdist(nu, beta, sJ3, sL1)
    alpha32 <- expdist(nu, beta, sJ3, sL2)
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

nubetaNM <- getnubeta(myorigarg$U, method="Nelder-Mead", control=list(trace=1, REPORT=1))
nubetaBFGS <- getnubeta(myorigarg$U, method="BFGS", control=list(trace=1, REPORT=1))

c(nubetaNM$value, nubetaBFGS$value)

cbind(nubetaNM$par, nubetaBFGS$par)

nubetaBFGS$par["nu"] / nubetaBFGS$par["beta"]


getnu <- function(alphatab, ...)
{
  
  f <- function(x)
  {
    nu <- x[1]
    sJ1 <- x[2]
    sJ2 <- x[3]
    sJ3 <- x[4]
    sL1 <- x[5]
    sL2 <- x[6]
    T11 <- alphatab[1,1] - expdist2(nu, sJ1, sL1)
    T12 <- alphatab[1,2] - expdist2(nu, sJ1, sL2)
    T21 <- alphatab[2,1] - expdist2(nu, sJ2, sL1)
    T22 <- alphatab[2,2] - expdist2(nu, sJ2, sL2)
    T31 <- alphatab[3,1] - expdist2(nu, sJ3, sL1)
    T32 <- alphatab[3,2] - expdist2(nu, sJ3, sL2)
    
    (T11^2+T12^2+T21^2+T22^2+T31^2+T32^2)/2
  }
  gradf <- function(x)
  {
    nu <- x[1]
    sJ1 <- x[2]
    sJ2 <- x[3]
    sJ3 <- x[4]
    sL1 <- x[5]
    sL2 <- x[6]
    alpha11 <- expdist2(nu, sJ1, sL1)
    alpha12 <- expdist2(nu, sJ1, sL2)
    alpha21 <- expdist2(nu, sJ2, sL1)
    alpha22 <- expdist2(nu, sJ2, sL2)
    alpha31 <- expdist2(nu, sJ3, sL1)
    alpha32 <- expdist2(nu, sJ3, sL2)
    T11 <- alphatab[1,1] - alpha11
    T12 <- alphatab[1,2] - alpha12
    T21 <- alphatab[2,1] - alpha21
    T22 <- alphatab[2,2] - alpha22
    T31 <- alphatab[3,1] - alpha31
    T32 <- alphatab[3,2] - alpha32
    
    M <- rbind(
      c(sL1-sJ1, sJ1-sL2, sJ2-sL1, sJ2-sL2, sJ3-sL1, sJ3-sL2), #w.r.t. nu
      c(-nu, -nu, 0, 0, 0, 0), #w.r.t. sJ1
      c(0, 0, -nu, -nu, 0, 0), #w.r.t. sJ2
      c(0, 0, 0, 0, -nu, -nu), #w.r.t. sJ3
      c(nu, 0, nu, 0, nu, 0), #w.r.t. sL1
      c(0, nu, 0, nu, 0, nu) #w.r.t. sL2
    )
    M <- M %*% diag(c(T11, T12, T21, T22, T31, T32))
    M <- M %*% diag(c(alpha11, alpha12, alpha21, alpha22, alpha31, alpha32))
    M %*% rep(1, 6) #sum the six terms
  }
  x0 <- c("nu"=1, "sJ"=1:3/4, "sL"=7:8/8)
  cat("x0", x0, "\n")
  cat(f(x0), "\n")
  print(gradf(x0))
  
  matU <- rbind(
    c(1, rep(0, 5)),
    cbind(0, rep(-1, 2), rep(0, 2), rep(0, 2), diag(2)),
    cbind(0, rep(0, 2), rep(-1, 2), rep(0, 2), diag(2)),
    cbind(0, rep(0, 2), rep(0, 2), rep(-1, 2), diag(2))
  )
  print(matU)
  
  constrOptim(x0, f, grad=gradf, ui=matU, ci=rep(0, NROW(matU)), ...)
}

nuBFGS <- getnu(myorigarg$U, method="BFGS", control=list(trace=1, REPORT=1))



nubetaBFGS$par["nu"] / nubetaBFGS$par["beta"]
nuBFGS$par

c(nuBFGS$value, nubetaBFGS$value)


dStationJoueur1 <- function(x) 
  expdist(nubetaBFGS$par["nu"], nubetaBFGS$par["beta"], nubetaBFGS$par["sJ1"], x)
dStationJoueur2 <- function(x) 
  expdist(nubetaBFGS$par["nu"], nubetaBFGS$par["beta"], nubetaBFGS$par["sJ2"], x)

curve(dStationJoueur1(x), from = 0, to = 1, n =501, xlab=expression(s),
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
      xlab=expression(s), ylab=expression(bar(q)(s)), from=0, to=1, n=1001)
abline(v=nubetaBFGS$par[3:5], lty=2, col="grey")
abline(h=100, col="red")

curve(qbarf(x, nubetaNM$par, c(.5, .25, .75), resbase[1,]$par[1:3]),
      from=0, to=1, n=1001)
abline(v=nubetaNM$par[3:5], lty=2, col="grey")
abline(h=100, col="red")


