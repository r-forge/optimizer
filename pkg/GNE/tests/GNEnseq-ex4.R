
require(GNE)	



#-------------------------------------------------------------------------------
# (4) Example of GNE with 4 solutions(!)
#-------------------------------------------------------------------------------



F <- function(z, phi=phiMin, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	# cat("x", x, "\n")
	c(	2*(x[1] - 2)*(x[2]-4)^4 + lambda[1],
		2*(x[2] - 3)*x[1]^4 + lambda[2],
		phi(lambda[1], 1-sum(x), ...),
		phi(lambda[2], 2-2*x[1]-x[2], ...)
	)
}


JacF <- function(z, gphia, gphib, ...)
{
	x <- z[1:2]
	lambda <- z[3:4]
	idga1 <- gphia(lambda[1], 1- sum(x), ...)
	idgb1 <- gphib(lambda[1], 1- sum(x), ...)
	idga2 <- gphia(lambda[2], 2-2*x[1]-x[2], ...)
	idgb2 <- gphib(lambda[2], 2-2*x[1]-x[2], ...)

rbind(
c(2*(x[2]-4)^4, 8*(x[1] - 2)*(x[2]-4)^3, 1, 0),
c(8*(x[2] - 3)*x[1]^3, 2*x[1]^4, 0, 1),
c(-idgb1, -idgb1, idga1, 0),
c(-2*idgb2, -idgb2, 0, idga2)
)	
}


#list of true GNEs
trueGNE <- rbind(c(2, -2, 0, 5*2^5),
	c(-2, 3, 8, 0),
	c(0, 1, 4*3^4, 0),
	c(1, 0, 2^9, 6))
colnames(trueGNE) <- c("x1", "x2", "lam1", "lam2")
rownames(trueGNE) <- 1:4

#check GNE
apply(trueGNE, 1, function(x) sqrt(sum(F(x)^2)) )

#is the jacobian singular?
for(i in 1:4)
{ 
	cat(">>>\t", i, "\n")
	print(JacF(trueGNE[i, ], gphia= GrAphiMin, gphib= GrBphiMin))
print( try( solve(JacF(trueGNE[i, ], gphia= GrAphiMin, gphib= GrBphiMin)) ) )
	cat("\n")
}


#check constraint
apply(trueGNE[, 1:2], 1, function(x) 1-sum(x) )
apply(trueGNE[, 1:2], 1, function(x) 2-2*x[1]-x[2] )


plot(trueGNE[, "x1"], trueGNE[, "x2"], pch=4, col="red")
grid()


#a simple test
z0 <- c(10, 10, 1, 1)

GNE.nseq(z0, F, JacF, list(phi=phiMin), list(gphia= GrAphiMin, gphib= GrBphiMin), method="Newton")

GNE.nseq(z0, F, JacF, list(phi= phiFB), list(gphia= GrAphiFB, gphib= GrBphiFB), method="Newton")

#-------------------------------------------------------------------------------
#random initial points
n <- 20
set.seed(1234)
initpt <- cbind(runif(n, -10, 10), runif(n, -10, 10), 1, 1)

NewLnsrch <- function(i, echo=FALSE)
{ 
	if(echo)
		cat("______", initpt[i, ], "\n")
	res <- GNE.nseq(initpt[i, ], F, JacF, list(phi=phiMin), list(gphia= GrAphiMin, gphib= GrBphiMin), method="Newton") 
	if(echo)
		print(res$par)
	c(res$par, res$value )
}

respt <- t( sapply(1:NROW(initpt), NewLnsrch ) )

totalres <- cbind(1:NROW(initpt), initpt[, 1:2], NA, round(respt, 3))

#remove non convergent optimization sequences
finalres <- totalres[totalres[,9] == 0, -9]
colnames(finalres) <- c("num", "x1 init", "x2 init", "", "x1 final", "x2 final", "lam 1", "lam 2")

print(finalres)
print(trueGNE)

?GNE


#-------------------------------------------------------------------------------
#test LM
z0 <- c(10, 10, 1, 1)
z0 <- c(-4, 4, 1, 1)

#min function
resMINLM <- bench.GNE.nseq.LM(z0, F, JacF, argPhi=list(phi=phiMin), argjac=list(gphia= GrAphiMin, gphib= GrBphiMin), control=list(maxit=300, ftol=1e-7))

resMINLM$compres

resMINLM$reslist

#FB function
resFBLM <- bench.GNE.nseq.LM(z0, F, JacF, argPhi=list(phi= phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), control=list(maxit=300, ftol=1e-7))

resFBLM$compres

GNE.nseq(z0, F, JacF, argPhi=list(phi= phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), method="Levenberg-Marquardt", global="gline", control=list(LM.param="min", trace=1))

GNE.nseq(z0, F, JacF, argPhi=list(phi= phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), method="Levenberg-Marquardt", global="qline", control=list(LM.param="min", trace=3))


GNE.nseq(z0, F, JacF, argPhi=list(phi= phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), method="Levenberg-Marquardt", control=list(LM.param="adaptive", trace=1, maxit=1000))



plot(trueGNE[, "x1"], trueGNE[, "x2"], pch=4, col="red", xlim=c(-2, 10), ylim=c(-2, 10))
points(z0[1], z0[2], pch=3)
grid()



#-------------------------------------------------------------------------------
#benchmark
z0 <- c(10, 10, 1, 1)
# z0 <- c(1/2, 1/2, 0, 0)

# F(z0)
# JacF(z0, gphia= GrAphiMin, gphib= GrBphiMin)

#min function
resMin <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMin), argjac=list(gphia= GrAphiMin, gphib= GrBphiMin), echo=FALSE)

resMin$compres
resMin$reslist[[1]]



#FB function
resFB <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiFB), argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), echo=FALSE)

resFB$compres



#Mangasarian function
resMan <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMan, f=function(t) t^3), argjac=list(gphia= GrAphiMan, gphib= GrBphiMan, fprime=function(t) 3*t^2), echo=FALSE, control=list(maxit=200))

resMan$compres

#LT function
resLT <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiLT, q=4), argjac=list(gphia= GrAphiLT, gphib= GrBphiLT, q=4))

resLT$compres


#KK function
resKK <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiKK, lambda=3/2), argjac=list(gphia= GrAphiKK, gphib= GrBphiKK, lambda=3/2))

resKK$compres

	