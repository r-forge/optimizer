require(optimz)
source("~/grose.R")
xx <- rep(pi,4)
ans <- Rvmmin(xx, grose.f)
ans
print(ans)
anso <- optimx(xx, grose.f, method="Rvmmin")
print(anso)