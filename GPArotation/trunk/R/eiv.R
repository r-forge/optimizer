eiv <- function(L,  identity=seq(ncol(L)), normalize=FALSE){
   A1 <- L[ identity, , drop=FALSE]
   Th <- solve(A1)
   B <- array(NA, dim(L))
   B[identity, ] <- diag(1, nrow(identity))
   B[-identity,] <- L[-identity,, drop=FALSE] %*% Th
   list(loadings=B, Th=Th, method="eiv", orthogonal=FALSE, convergence=TRUE,
        Phi= t(Th) %*% Th)
   }

