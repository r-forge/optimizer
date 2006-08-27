eiv <- function(L,  identity=seq(ncol(L)), normalize=FALSE){
   A1 <- L[ identity, , drop=FALSE]
   Th <- solve(A1)
   if(1e-14 < max(abs(diag(1, length(identity)) - A1 %*% Th)))
      warning("Inverse is not well conditioned. Consider setting identity to select different rows.")
   B <- array(NA, dim(L))
   B[identity, ] <- diag(1, length(identity))
   B[-identity,] <- L[-identity,, drop=FALSE] %*% Th
   dimnames(B) <- list(dimnames(L)[[1]], paste("factor", seq(length(identity))))
   list(loadings=B, Th=Th, method="eiv", orthogonal=FALSE, convergence=TRUE,
        Phi= t(Th) %*% Th)
   }

