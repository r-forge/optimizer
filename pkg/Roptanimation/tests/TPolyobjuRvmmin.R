## @knitr polyexample

# library(Roptanimation)

# source("../R/smallpoly.R")

# Example code -- seems to work for nv=6, but not otherwise ??
# nvex <- as.numeric(readline("Number of vertices ="))
nvex = 6 # For the automated R CMD check
# Note the as.numeric
nv <- nvex # to copy value above
nsave <- 2*nv-1
cat("There are ",nv," vertices\n")

cat("Polygon data:\n")
reghex <- polysetup(nv, defsize=1)
regxy <- polypar2XY(reghex$par0)
regarea <- polyarea(reghex$par0)
cat("reghex area =", regarea, "\n")
myhex <- polysetup(nv, defsize=0.98)
# Note: default "size" is 0.98, not 1
print(myhex)
cat("Area:\n")
myhexa <- polyarea(myhex$par0)
print(myhexa)
cat("XY coordinates\n")
myheXY <- polypar2XY(myhex$par0)
print(myheXY)
# plot(myheXY$x, myheXY$y, type="l")
cat("Constraints:\n")
myhexc<-polydistXY(myheXY)
print(myhexc)
cat("Vertex distances:")
print(sqrt(myhexc))

start <- myhex$par0
# options(scipen=5) # did not work

pt1 <- new.env()
pt1$psave <- matrix(c(reghex$par0, -1, regarea), nrow=1)
pt1$besta <- 0 # best area
pt1$regarea <- regarea
#- ?? names!
f0 <- polyobju(start)


library(optimrx)
# ub <- c(rep(1,(nv-1)), rep(0.75*pi, (nv-2))) # approx for angles
# lb <- c(rep(0, (2*nv-3)))
# sol <- optimr(start, polyobju, polygradu, method="Rvmmin", lower=lb, upper=ub,

sol <- optimr(start, polyobju, polygradu, method="Rvmmin",
                control=list(trace=1, maxit=1000), penfactor=1e-5)
sol

