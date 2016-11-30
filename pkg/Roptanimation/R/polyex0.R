## @knitr polyex0

# Example code
nv <- 6
cat("Polygon data:\n")
myhex <- polysetup(nv)
print(myhex)
x0 <- myhex$par0 # initial parameters
cat("Area:\n")
myhexa <- polyarea(x0)
print(myhexa)
cat("XY coordinates\n")
myheXY <- polypar2XY(x0)
print(myheXY)
plot(myheXY$x, myheXY$y, type="l")
cat("Constraints:\n")
ub <- c(rep(1,(nv-1)), rep(0.75*pi, (nv-2))) # approx for angles
lb <- c(rep(0, (2*nv-3)))
myhexc<-polydistXY(myheXY)
print(myhexc)
cat("Vertex distances:")
print(sqrt(myhexc))
cat("check distances with polypar2distXY\n")
try1 <- polypar2distXY(x0)
print(try1)
cat("check distances with polypardist2 augmenting output with parameter squares\n")
try2 <- polypardist2(x0)
try2 <- c(x0[1:(nv-1)]^2, try2)
print(try2)
cat("Max abs difference = ",max(abs(try1-try2)),"\n")

