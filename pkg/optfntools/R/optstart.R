# optstart.R
optstart<-function(npar){ # create structure for optimization computations
   # npar is number of parameters
   if (is.null(npar)) stop("optstart: You must specify npar")
   if (npar<=0) stop("optstart: npar must be a positive integer")
   if (npar != round(npar)) stop("optstart: npar must be a positive integer")
   opdata<-list(MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
       KFN=0, KGR=0, KHESS=0)
   # may be other stuff -- bounds, masks??
   localOPCON<-list2env(opdata)
}

