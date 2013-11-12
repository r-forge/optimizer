ztime<-function(x,ipivot) {
   ## ---------------------------------------------------------
   ##  this routine multiplies the vector x by the 
   ##  constraint matrix z
   ## ---------------------------------------------------------
   ind <- which(ipivot!=0);
   x[ind] <- 0
   x1 <- x;
}