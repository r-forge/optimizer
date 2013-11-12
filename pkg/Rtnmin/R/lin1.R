lin1 <- function(p, x, f, alpha, g, sfun, ...){
   ## ---------------------------------------------------------
   ##  line search (naive)
   ## ---------------------------------------------------------
   ##  set up
   ## ---------------------------------------------------------
##   cat("lin1: alpha=",alpha,"  p:\n")
##   print(p)

   if (is.null(alpha)) alpha <- 0

   ierror <- 3 
   xnew   <- x 
   fnew   <- f 
   gnew   <- g 
   maxit  <- 15 
   if (alpha == 0) {
      ierror <- 0
      maxit <- 1
   }
   alpha1 <- alpha 
   ## ---------------------------------------------------------
   ##  line search
   ## ---------------------------------------------------------
   for (itcnt in 1:maxit) {
      xt <- x + alpha1*p 
      fg <- sfun(xt)
      ft<-fg$f
      gt<-fg$g # may simplify later
      if (ft < f) {
         ierror <- 0
         xnew   <- xt
         fnew   <- ft
         gnew   <- gt
         ## cat("about to break in lin1\n")
         break
      }
      alpha1 <- alpha1 / 2
   }
   if (ierror == 3) { alpha1 <- 0 }
##   cat("lin1: itcnt=",itcnt,"\n")
   nf1 <- itcnt 

   ## never used in SGN code
   ## if (nargout == 7) { ## ?? what is nargout?
   ##  dfdp <- as.numeric(crossprod(gt, p) )
   ##  varargout{1} <- dfdp ## ?? what is varargout
   ## end
   result<-list(xnew=xnew, fnew=fnew, gnew=gnew, nf1=nf1,
          ierror=ierror, alpha1=alpha1)
}
