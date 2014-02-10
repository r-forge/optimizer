library(GNE)

#-------------------------------------------------------------------------------
# (2) Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------



#constants
mysmallarg <- list(d= 20, lambda= 4, rho= 1)
dimx <- c(1, 1)


#O_i(x)
obj <- function(x, i, arg)
{
	-(arg$d - arg$lambda - arg$rho*(x[1]+x[2]))*x[i]
}
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	res <- -arg$rho * x[i]
	if(i == j)
		res <- res + arg$d - arg$lambda - arg$rho*(x[1]+x[2])
	-res
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
	arg$rho * (i == j) + arg$rho * (j == k)	


dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	-x[i]
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	-1*(i == j)
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0

dimmu <- 1	
#joint function
h <- function(x)
	-min(x[1], x[2])	
grh <- function(x, j)
	-c(1*(x[1] < x[2]), 1*(x[2] < x[1]))[j]
heh <- function(x, j, k)
	0	



z0 <- rep(0, sum(dimx)+sum(dimlam))

#true value is (16/3, 16/3, 0, 0) 
myGNE <- GNE.nseq(z0, dimx, dimlam, grobj=grobj, mysmallarg,
	heobj=heobj, mysmallarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton")
	


#true value is (8, 4) 
mySE1 <- SE.nseq(1, c(2/3, 2/3, 0, 0), dimx, dimlam, 
	obj=obj, argobj=mysmallarg,
	grobj=grobj, arggrobj=mysmallarg,
	heobj=heobj, argheobj=mysmallarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, 
	control.leader=list(trace=TRUE), silent=TRUE, global="dbldog")

#true value is (8, 4) 
mySE2 <- SE.nseq(1, c(2/3, 2/3, 0), dimx,  
	obj=obj, argobj=mysmallarg,
	grobj=grobj, arggrobj=mysmallarg,
	heobj=heobj, argheobj=mysmallarg, 
	joint=h, grjoint=grh, hejoint=heh, dimmu=1, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, 
	control.leader=list(trace=TRUE), silent=FALSE, global="dbldog")


#-------------------------------------------------------------------------------
#comparison between NE and SE
	
res <- cbind(GNE=myGNE$par[1:2], objGNE=
	-c(obj(myGNE$par[1:2], 1, mysmallarg), obj(myGNE$par[1:2], 2, mysmallarg)),
	SE=mySE1$par, objSE=-c(obj(mySE1$par, 1, mysmallarg), obj(mySE1$par, 2, mysmallarg)))
rownames(res) <- c("P1","P2")
print(res)	
	
	
	
#-------------------------------------------------------------------------------
#manual computation
	
gbis <- function(xm1, i, arg)
	g(c(arg$x1, xm1),i+1)
grgbis <- function(xm1, i, j, arg)
	grg(c(arg$x1, xm1), i+1, j+1)
hegbis <- function(xm1, i, j, k, arg)
	heg(c(arg$x1, xm1), i+1, j+1, k+1)
		
grobjbis <- function(xm1, i, j, arg)
	grobj(c(arg$x1, xm1), i+1, j+1, arg$add)
heobjbis <- function(xm1, i, j, k, arg)
	heobj(c(arg$x1, xm1), i+1, j+1, k+1, arg$add)

bestresponse <- function(x1)
{
	myarg <- list(x1=x1, add=mysmallarg)
	GNE.nseq(rep(1, 2), dimx[-1], dimlam[-1], grobjbis, arggrobj=myarg, 
		heobjbis, argheobj=myarg, gbis, argconstr=myarg, grgbis,
		arggrconstr=myarg, hegbis, argheconstr=myarg,
		compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
		control=list(trace=0))$par[1]
}

bestresponse(2/3) #should be 7+2/3

#to be minimized
objleader <- function(x, arg)
{
	bestresp <- bestresponse(x)
	-(arg$d - arg$lambda - arg$rho*(x+bestresp))*x
}


res2 <- optim(0, objleader, arg=list(d= 20, lambda= 4, rho= 1), method="BFGS", control=list(trace=1))

#true value (8,4)
xstar <- c(res2$par, bestresponse(res2$par))
xstar

#-------------------------------------------------------------------------------
#graphical check

xlead <- seq(1, 15, length=101)
fnofxlead <- SE.objleaders(xlead, 1, c(2/3, 2/3, 0, 0), dimx, dimlam, 
	obj=obj, argobj=mysmallarg,
	grobj=grobj, arggrobj=mysmallarg,
	heobj=heobj, argheobj=mysmallarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, 
	silent=FALSE)

plot(xlead, fnofxlead, type="l")
abline(v=mySE1$par[1], col="red")
