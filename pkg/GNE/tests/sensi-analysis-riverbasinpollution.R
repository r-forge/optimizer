if(!require("GNE"))stop("this test requires package GNE.")

#if(FALSE)
{
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
  
  dimlam <- c(2, 2, 2)
  #g_i(x)
  g <- function(x, i, arg)
    c(sum(arg$U[, 1] * arg$E * x[1:3]) - arg$K[1],
      sum(arg$U[, 2] * arg$E * x[1:3]) - arg$K[2],
      -x[1],
      -x[2],
      -x[3])
  #Gr_x_j g_i(x)
  grg <- function(x, i, j, arg)
    c(arg$U[j, 1] * arg$E[j], 
      arg$U[j, 2] * arg$E[j], 
      -1*(i ==j), 
      -1*(i ==j), 
      -1*(i ==j))
  #Gr_x_k Gr_x_j g_i(x)
  heg <- function(x, i, j, k, arg)
    c(0, 0, 0, 0, 0)
  
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
  
  
  getNE <- function(x, control=list(maxit=100, trace=0), check=TRUE, output="equilibrium")
  {
    output <- match.arg(output, c("equilibrium", "fullobject"))
    res <- sapply(1:NROW(x), function(i)
    {
      myarg <- list(
        C = cbind(x[i,paste("C",1:3,sep="")], x[i,paste("C",1:3+3,sep="")]),
        U = cbind(x[i,paste("U",1:3,sep="")], x[i,paste("U",1:3+3,sep="")]),
        K = x[i,paste("K",1:2,sep="")],
        E = x[i,paste("E",1:3,sep="")],
        D = x[i,paste("D",1:2,sep="")]
      )
      z0 <- rep(1, sum(dimx)+sum(dimlam))
      
      res <- GNE.nseq(z0, dimx, dimlam, grobj=grobj, arggrobj=myarg, heobj=heobj, argheobj=myarg, 
                      constr=g, myarg, grconstr=grg, myarg, heconstr=heg, myarg, 
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
  
  getNE2 <- function(x, control=list(maxit=100, trace=3), check=TRUE, output="equilibrium")
  {
    output <- match.arg(output, c("equilibrium", "fullobject"))
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
                      control=control, silent=FALSE)
      
      if(output == "equilibrium")
        return(res$par[1:3])
      else
        return(res)
    }
    )
    t(res)
  }
  
  
  n <- 100
  nbparam <- 19
  Xinputmatrix <- t(replicate(n, unlist(myarg0) * (1 + rnorm(nbparam, 0, .1))))
  X2inputmatrix <- t(replicate(n, unlist(myarg0) * (1 + rnorm(nbparam, 0, .1))))
  
  #pairs(Xinputmatrix[, 1:10])
  
  basesetting <- matrix(unlist(myarg0), nrow=1, dimnames = list(NULL, names(unlist(myarg0))))
  
  resbase <- getNE(basesetting, output = "full")[1,]
  resbase2 <- getNE2(basesetting, output= "full")[1,]
  
  cbind(resbase$par[1:3], resbase2$par[1:3])
  cbind(resbase$value, resbase2$value)
  cbind(resbase$fvec[1:3], resbase2$fvec[1:3])
  c(resbase$iter, resbase2$iter)
  cbind(resbase$counts, resbase2$counts)
  
  Youtputs <- t(getNE(Xinputmatrix))
  idx <- ridx <- ridx <- rowSums(is.na(Youtputs)) == 0
  Youtputs <- Youtputs[idx, ]
  Xinputmatrix <- Xinputmatrix[idx, ]
  X2inputmatrix <- X2inputmatrix[idx, ]
  
  Yarithmean <- apply(Youtputs, 1, mean)
  Ygeomean <- apply(Youtputs, 1, function(x) prod(x^(1/3)))
  
  # plot(Yarithmean)
  # points(Ygeomean, col="red")
  # hist(Ygeomean)
  # hist(Yarithmean)
  
  
  
  # packages
  library(sensitivity)
  
  respccArith <- pcc(data.frame(Xinputmatrix), Yarithmean)
  respccArith$PCCC
  respccGeo <- pcc(data.frame(Xinputmatrix), Ygeomean)
  respccGeo$PCCC
  
  
  ressrcArith <- src(data.frame(Xinputmatrix), Yarithmean)
  ressrcGeo <- src(data.frame(Xinputmatrix), Ygeomean)
  
  # resobol <- sobol(getNE, Xinputmatrix, X2inputmatrix, check=FALSE)
  
}