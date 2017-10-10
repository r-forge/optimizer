library(optimrx)
sall <-  opm(x0, polyobj, polygrad, method="ALL", lower=lb, upper=ub, 
             control=list(trace=1, kkt=FALSE), penfactor=1e-5)
sall <- summary(sall, order=value)
print(sall)
best1 <- coef(sall)[1,]
polyarea(best1)