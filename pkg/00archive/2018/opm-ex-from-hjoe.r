# logistic regression, different optimizers in optimrx

# minimize negative log-likelihood = nllk
# y is an nx1 binary 0-1 vector
# xdat is a nxp data matrix
# param is a (p+1) dimensional vector, including an intercept
logregnllk=function(param,xdat,y)
{ b0=param[1]
  bvec=param[-1]
  if(is.vector(xdat)) { tem=b0+bvec[1]*xdat }
  else { tem=b0+xdat%*%bvec }
  nllk= -sum(tem*y) + sum(log(1+exp(tem)))
  nllk
}

library(MASS)
data(Pima.tr)

fit2l=glm(type~glu+bmi+ped+age,family=binomial,data=Pima.tr)
print(summary(fit2l))
#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -9.971388   1.527587  -6.528 6.69e-11 ***
#glu          0.031255   0.006627   4.716 2.40e-06 ***
#bmi          0.077030   0.032251   2.388 0.016921 *  
#ped          1.719794   0.656088   2.621 0.008760 ** 
#age          0.058603   0.017574   3.335 0.000854 ***
#    Null deviance: 256.41  on 199  degrees of freedom
#Residual deviance: 181.08  on 195  degrees of freedom
#AIC: 191.08

xdat=as.matrix(Pima.tr[c('glu','bmi','ped','age')])
xcenter=scale(xdat,center=T,scale=F)
y=as.numeric(Pima.tr$type)-1
sub=data.frame(cbind(xcenter,y))

# centered covariates
fit2c=glm(y~glu+bmi+ped+age,family=binomial,data=sub)
print(summary(fit2c))
#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -0.933686   0.195641  -4.772 1.82e-06 ***
#glu          0.031255   0.006627   4.716 2.40e-06 ***
#bmi          0.077030   0.032251   2.388 0.016921 *  
#ped          1.719794   0.656088   2.621 0.008760 ** 
#age          0.058603   0.017574   3.335 0.000854 ***
#    Null deviance: 256.41  on 199  degrees of freedom
#Residual deviance: 181.08  on 195  degrees of freedom
#AIC: 191.08

ybar=mean(y)
b0=log(ybar/(1-ybar))

library(optimrx)
cat("\n\nopm try 1\n")
mleall=opm(c(b0,0,0,0,0),logregnllk, method="ALL",hessian=T,xdat=xdat,y=y,
  control=list(kkt=F))
summary(mleall, order=value)
# nlm and other fails, but nlm works find when used on its own

cat("\n\nopm with better starting point\n")
mleallc=opm(c(b0,0,0,0,0),logregnllk,  method="ALL",hessian=T,xdat=xcenter,y=y,
  control=list(kkt=F))
summary(mleallc, order=value)
# still a surprising number of optimizers fail for this better start

mleallnd=optimr(c(b0,0,0,0,0),logregnllk,gr="grnd", method="nlm",control=list(trace=1),xdat=xcenter,y=y)
print(mleallnd)

mleallnc=optimr(c(b0,0,0,0,0),logregnllk,gr="grcentral", method="nlm",control=list(trace=1),xdat=xcenter,y=y)
print(mleallnd)

mleallnf=optimr(c(b0,0,0,0,0),logregnllk,gr="grfwd", method="nlm",control=list(trace=1),xdat=xcenter,y=y)
print(mleallnd)

mleallnlm=nlm(logregnllk, c(b0,0,0,0,0), xdat=xcenter, y=y, print.level=1)
print(mleallnlm)
