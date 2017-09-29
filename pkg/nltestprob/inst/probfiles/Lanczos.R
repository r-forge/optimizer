# Optimization test function lanczos1
# lanczos1 from NISTnls
# ??ref...

# Only from res, so remove
# lanczos1.f <- function(x) {
#   res<-lanczos1.res(x)
#   f<-sum(res*res)
# }

lanczos.res <- function(b, pdata) {
   xx<-pdata$x # Assume we have properly structured dataframe in pdata
   yy<-pdata$y
   res <- rep(NA, length(xx)) # just in case
   b1<-b[1] # get parameters from the parameter vector
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   res<-as.vector(b1*exp(-b2*xx) + b3*exp(-b4*xx) + b5*exp(-b6*xx) - yy)
}

# lanczos - Jacobian
lanczos.jac <- function(b, pdata) {
   # 170303
   xx<-pdata$x
   yy<-pdata$y
   n<-length(b)
   m<-length(xx)
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   b4<-b[4]
   b5<-b[5]
   b6<-b[6]
   J<-matrix(0,m,n) # define the size of the Jacobian and zero it
   J[,1] <- exp(-b2*xx)
   J[,3] <- exp(-b4*xx)
   J[,5] <- exp(-b6*xx)
   J[,2] <- -b1*xx*exp(-b2*xx)
   J[,4] <- -b3*xx*exp(-b4*xx)
   J[,6] <- -b5*xx*exp(-b6*xx)
   attr(J,"gradient") <- J
   J
}

lanczos1.h <- function(x) {
stop("not defined")
   JJ<-lanczos1.jac(x)
   H <- t(JJ) %*% JJ
   res<-lanczos1.res(x)
stop("not defined")

}

lanczos1.g<-function(x) {
#   stop("not defined")
   JJ<-lanczos1.jac(x)
   res<-lanczos1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

lanczos1.fgh<-function(x) {
   f<-lanczos1.f(x)
   g<-lanczos1.g(x)
   H<-lanczos1.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

probdata<-function(dataname, pkgname) {
   library(pkgname, character.only=TRUE) # get parent collection
   eval(parse(text=data(list=dataname))) # and load up the data into x and y
}


lanczos.setup <- function(idxdata=NULL, idxstart=NULL) {
   if(is.null(idxdata)) {
      cat("There are three Lanczos data sets, indexed 1:3\n")     
   } else {
      if (idxdata %in% 1:3) {
          dsetname <- paste("Lanczos",idxdata,sep='')
          mypdata <- probdata(dsetname, "NISTnls")
      } else { stop("ERROR-Lanczos data set index out of range = ",idxdata)}
   }
  if(is.null(idxstart)) {
    cat("There are two Lanczos starting points, indexed 1:2, or 0 to enter from keyboard\n")     
  } else {
    if (idxstart %in% 0:2) {
        if (idxstart == 0) {
             cat("Enter the 6 elements of the starting vector:")
             for (jj in 1:6) {
                  start[jj] <- readline("? ")
             }
        } else { 
             if (idxstart == 1) { 
                 start<-c(1.2,0.3,5.6,5.5,6.5,7.6)
             } else { 
                 start<-c(0.5,0.7,3.6,4.2,4,6.3) 
                 }
        }
    } else { stop("ERROR-Lanczos starting point index out of range = ",idxdata)}
  }
  Lz <- list(mypdata=mypdata, start=start)
}
