## @knitr ##Petran.prb
# This is file ##Petran.prb
rm(list=ls())
probname <- "##Petran"
probdesc <- "

  File is
  pH	mumax
  4	0
  4.1	0
  4.3	0
  4.5	0
  4.7	0.112
  5	0.229
  6	0.8
  7	0.93
  8	0.83
  9	0.285
  9.2	0.232
  9.4	0
  9.6	0
  9.8	0
  10	0

"

#- Note: environment / list "counters" must already exist

if (exists("pe")) { 
  rm("pe")  
}

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

mypdata <- read.csv("petran.csv")

#- nls format expression
ratkowsky<-as.formula("mumax ~ b*( (pH-pHmin) * (1-exp(c*(pH-pHmax))) )^2")

#- setup
### Petran Ratkowsky model
### describing a bacterial specific growth rate as
### a function of pH

# fit with nls
test<-try(fit<-nls(ratkowsky, mypdata,list(b=10,c=0.01,pHmin=4,pHmax=10)))
fit
test<-try(fit2<-nls(ratkowsky, mypdata,list(b=10,c=0.01,pHmin=4,pHmax=10),algorithm="port"))
fit2
test<-try(fit3<-nls(ratkowsky, mypdata,list(b=16,c=0.03,pHmin=4,pHmax=10),algorithm="port"))
fit3
x11()
plotfit(fit,smooth=TRUE)
# confint(fit) # ne tourne pas
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
	mumaxtheo<-b*((mypdata$pH-pHmin)*(1-exp(c*(mypdata$pH-pHmax))))^2
	RSS<-sum((mypdata$mumax-mumaxtheo)^2)
	return(RSS)
}
res1<-optim(par=c(10,0.01,4,10),fn=RSS) # Nelder and Mead 
res1
x11()
preview(ratkowsky, mypdata,list(b=13.97,c=-0.031,pHmin=4.2,pHmax=9.9))
res2<-optim(par=c(10,0.01,4,10),fn=RSS,method="BFGS") # quasi-Newton 
res2
res3<-optim(par=c(10,0.01,4,10),fn=RSS,method="CG") # conjugate gradients  
res3
res4<-optim(par=c(10,0.01,4,10),fn=RSS,method="SANN") # variant of simulated annealing  
res4
