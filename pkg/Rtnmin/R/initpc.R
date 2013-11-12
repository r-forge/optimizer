initpc <- function(d, upd1, ireset) {
##---------------------------------------------------------
## initialize the diagonal preconditioner (d -- a vector!)
## ---------------------------------------------------------
## global hyk sk yk sr yr yksk yrsr ?? do we have these
#  global vectors hyk sk yk sr yr & scalars yksk yrsr
## ---------------------------------------------------------
#%   cat("In initpc -- sk and hyk if upd1 false: ", upd1,"\n")
#%   print(envjn$sk)
   if (upd1) { 
      td <- d
   } else {
      if (ireset) {
         envjn$hyk  <- d * envjn$sk # vector * vector by element??
         sds  <- as.numeric(crossprod(envjn$sk, envjn$hyk)) # matrix multiply
         if (all(envjn$hyk == 0)) { cat("INITPC: envjn$hyk = 0 \n") }
         if (sds == 0) { cat("INITPC: sds = 0 \n") }
         td   <- d - d*d*envjn$sk*envjn$sk/sds + envjn$yk*envjn$yk/envjn$yksk 
         # by element 
      } else {
         envjn$hyk  <- d * envjn$sr # vector * vector
         sds  <- as.numeric(crossprod(envjn$sr, envjn$hyk))
         srds <- as.numeric(crossprod(envjn$sk, envjn$hyk))
         yrsk <- as.numeric(crossprod(envjn$yr, envjn$sk))
         envjn$hyk  <- d*envjn$sk - envjn$hyk*srds/sds + envjn$yr*yrsk/envjn$yrsr
         td   <- d - d*d*envjn$sr*envjn$sr/sds+envjn$yr*envjn$yr/envjn$yrsr
         sds  <- as.numeric(crossprod(envjn$sk, envjn$hyk))
         td   <- td - envjn$hyk*envjn$hyk/sds + envjn$yk*envjn$yk/envjn$yksk
      }
   }
#%   cat("td:")
#%   print(td)
   ans<-list(td=td)
}
