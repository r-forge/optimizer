## scalsqr.R -- a scaled quadratic function

rawfn<-function(par) { # unscaled function
	npar<-length(par)
	fval<-NA
	sq<-1:npar
	zz<-par-sq
	if (npar == 1) {
		fval<-zz*zz
	} else {
		fval<-0
		for (i in 1:(npar-1)) {
			fval<-fval+(zz[i+1]+zz[i])^2
		}
		if (npar>2) fval<-fval+(zz[1]+zz[npar])^2
	}
	return(fval)
}

scfn<-function(par) { # scaled function
	npar<-length(par)
	fval<-NA
	sq<-1:npar
	zz<-(par-sq)*(2^sq)
	if (npar == 1) {
		fval<-zz*zz
	} else {
		fval<-0
		for (i in 1:(npar-1)) {
			fval<-fval+(zz[i+1]+zz[i])^2
		}
		if (npar>2) fval<-fval+(zz[1]+zz[npar])^2
	}
	return(fval)
}


test1<-function() {
#   source("scalsqr.R")
   set.seed(123)
   npar<-readline("Number of parameters=")
   cat("rawfn\n")
   for (i in 1:25) {
       ss<-runif(npar,-20,20)
       au<-optim(ss,rawfn)
       for (j in 1:npar) {cat(au$par[j])}
       cat(" f=",au$value,"\n")
   }



}
