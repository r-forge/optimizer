library(nlstools)
d<-read.table("petran.txt",header=T)

# File is
# pH	mumax
# 4	0
# 4.1	0
# 4.3	0
# 4.5	0
# 4.7	0.112
# 5	0.229
# 6	0.8
# 7	0.93
# 8	0.83
# 9	0.285
# 9.2	0.232
# 9.4	0
# 9.6	0
# 9.8	0
# 10	0

### Petran Ratkowsky model
### describing a bacterial specific growth rate as
### a function of pH
ratkowsky<-as.formula("mumax ~ b*( (pH-pHmin) * (1-exp(c*(pH-pHmax))) )^2")

# fit with nls
test<-try(fit<-nls(ratkowsky,d,list(b=10,c=0.01,pHmin=4,pHmax=10)))
fit
test<-try(fit2<-nls(ratkowsky,d,list(b=10,c=0.01,pHmin=4,pHmax=10),algorithm="port"))
test<-try(fit3<-nls(ratkowsky,d,list(b=16,c=0.03,pHmin=4,pHmax=10),algorithm="port"))
fit3
x11()
plotfit(fit,smooth=TRUE)
confint(fit) # ne tourne pas
overview(fit)
cont<-nlsContourRSS(fit)
x11()
plot(cont,nlev=10,col=FALSE)


# fit using optim to minimize the RSS
RSS<-function(vpar)
{
	b<-vpar[1]
	c<-vpar[2]
	pHmin<-vpar[3]
	pHmax<-vpar[4]

	# calcul de la SCE
	mumaxtheo<-b*((d$pH-pHmin)*(1-exp(c*(d$pH-pHmax))))^2
	RSS<-sum((d$mumax-mumaxtheo)^2)
	return(RSS)
}
res1<-optim(par=c(10,0.01,4,10),fn=RSS) # Nelder and Mead 
res1
x11()
preview(ratkowsky,d,list(b=13.97,c=-0.031,pHmin=4.2,pHmax=9.9))
res2<-optim(par=c(10,0.01,4,10),fn=RSS,method="BFGS") # quasi-Newton 
res2
res3<-optim(par=c(10,0.01,4,10),fn=RSS,method="CG") # conjugate gradients  
res3
res4<-optim(par=c(10,0.01,4,10),fn=RSS,method="SANN") # variant of simulated annealing  
res4
