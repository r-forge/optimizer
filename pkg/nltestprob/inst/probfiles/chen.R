## Chen function ?? ref?

chen.f <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
   sum (res * res)
}

chen.g <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
   jj<-chen.jac(x)
   gg<- 2.0 * jj %*% as.vector(res)
   return(gg)
}

chen.res <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
}

chen.jac <- function(x) {
   n<-length(x)
   jj<-matrix(0.0, n, n)
   for (i in 1:n) {
     jj[i,i] <- 0.5 * (1.0/x[i] + exp(x[i])) * (1.0-v[i]/sqrt(v[i]^2 + 5e-04))
   } 
   return(jj)
}


chen.setup <- function(n=NULL, probcase=NULL) {
   if ( is.null(n) ) {
      n <- readline("Order of problem (n):")
   }
   if (is.null(probcase)) { probcase <- 1 }
   if (probcase == 1) { x<-rep(2,n) }
#  else if (probcase == 2) {       }
   lower<-rep(0.0, n) # Because log cannot take negative argument?
   upper<-rep(100.0, n)
   bdmsk<-rep(1,n)
   gsu<-list(x=x,lower=lower,upper=upper,bdmsk=bdmsk)
   return(gsu)
}






