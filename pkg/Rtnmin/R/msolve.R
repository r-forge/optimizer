msolve <- function(g, upd1, ireset, first, d) {
##---------------------------------------------------------
## This routine acts as a preconditioning step for the
## linear conjugate-gradient routine.  It is also the
## method of computing the search direction from the
## gradient for the non-linear conjugate-gradient code.
## It represents a two-step self-scaled bfgs formula.
##---------------------------------------------------------
#%    cat("msolve, envjn$john=",envjn$john,"\n")
   envjn$john<-"msolve"
   if (upd1) {
#%       cat("upd1 is TRUE\n")
      y <- g/d
   } else {
      gsk <- as.numeric(crossprod(g, envjn$sk))
      if (ireset) {
         envjn$hg <- g/d 
         if (first) {
	    envjn$hyk   <- envjn$yk/d 
            ykhyk <- as.numeric(crossprod(envjn$yk, envjn$hyk))
         }
         ghyk <- as.numeric(crossprod(g, envjn$hyk)) 
         y    <- ssbfgs(envjn$sk,envjn$hg,envjn$hyk,envjn$yksk,ykhyk,gsk,ghyk)
      } else {
         envjn$hg <- g/d 
         if (first) {
            envjn$hyk   <- envjn$yk/d 
            envjn$hyr   <- envjn$yr/d 
            envjn$yksr  <- as.numeric(crossprod(envjn$yk, envjn$sr))
            ykhyr <- as.numeric(crossprod(envjn$yk, envjn$hyr))
         }
         gsr <- as.numeric(crossprod(g, envjn$sr))
         ghyr <- as.numeric(crossprod(g, envjn$hyr))
         if (first) {
            yrhyr <- as.numeric(crossprod(envjn$yr, envjn$hyr))
         }
         envjn$hg <- ssbfgs(envjn$sr,envjn$hg,envjn$hyr,envjn$yrsr,yrhyr,gsr,ghyr)
         if (first) {
            envjn$hyk <- ssbfgs(envjn$sr,envjn$hyk,envjn$hyr,envjn$yrsr,yrhyr,envjn$yksr,ykhyr)
         }
         ykhyk <- as.numeric(crossprod(envjn$hyk, envjn$yk))
         ghyk  <- as.numeric(crossprod(envjn$hyk, g)) 
         y     <- ssbfgs(envjn$sk,envjn$hg,envjn$hyk,envjn$yksk,ykhyk,gsk,ghyk) 
      }
   }
   amsolve<-list(y=y) ## returns y
}

