## Test hessian computation in optimr
# thessian.R
library(optimr)

library(adagio)

x0 <- rep(pi,4)

thjn <- optimr(x0, fnRosenbrock, method="hjn", hessian=TRUE)
print(thjn)
thjn0 <-  optimr(x0, fnRosenbrock, method="hjn")
print(thjn0)

tlb <- optimr(x0, fnRosenbrock, grRosenbrock, method="L-BFGS-B", hessian=TRUE)
print(tlb)
tlb0 <-  optimr(x0, fnRosenbrock, grRosenbrock, method="L-BFGS-B")
print(tlb0)

tlbx <- optimr(x0, fnRosenbrock, method="L-BFGS-B", hessian=TRUE)
print(tlbx)
tlbx0 <-  optimr(x0, fnRosenbrock, method="L-BFGS-B")
print(tlbx0)

