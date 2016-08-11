require(BB)
set.seed(123) 

rosbkext.f <- function(x){
n <- length(x)
sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

n <- 20
p0 <- rnorm(n)

# 2 constraints: parameters sum to 1 and the last parameter is non-negative
# x[1] + ... + x[n] = 1
# x[n] >= 0
Amat <- rbind(rep(1,n), c(rep(0,n-1),1))
b <- c(1, 0)

#  with projectLinear as in release 2014.1-1 next gave
#  Failure in initial projection!Error in solve.QP(dvec = rep(0, n), Amat = t(A), 
#  bvec = b - c(A %*% par),  : 
#  Error in projection

ans <- spg(par=p0, fn=rosbkext.f, project="projectLinear", 
   projectArgs=list(A=Amat, b=b, meq=1), control=list(maxit=2500)) 


if(1e-3 < max(abs(ans$par -  1.0 ))){	
   print(ans$par, digits=18)
   cat("difference:\n")
   print(ans$par -  1.0 , digits=18)
   stop("converged to different parameter values!")
   }

if(1e-6 < max(abs(ans$value - 0.0 ))){
   print(ans$value, digits=18)
   stop("converged to different function value!")
   }

ans

# 2014 version and previous gave following but was indicating
#  convergence when it had not really been obtained.
# $par
#  [1]  5.460011e-01  3.001331e-01  9.300775e-02  1.186805e-02  3.388513e-03
#  [6]  3.259553e-03  3.258705e-03  3.258699e-03  3.258702e-03  3.258695e-03
# [11]  3.258700e-03  3.258702e-03  3.258694e-03  3.258701e-03  3.258699e-03
# [16]  3.258699e-03  3.258699e-03  3.258562e-03  3.237670e-03 -6.098637e-19
# 
# $value
# [1] 17.41526
# 
# $gradient
# [1] 0.0001785687
# 
# $fn.reduction
# [1] 5362.369
# 
# $iter
# [1] 52
# 
# $feval
# [1] 53
# 
# $convergence
# [1] 0
# 
# $message
# [1] "Successful convergence"
