ydat<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat<-1:length(ydat) # for testing
weeddata<-data.frame(y=ydat, t=tdat)
require(nlmrt)
require(minpack.lm)
hobmod <- "y~100*b1/(1+10*b2*exp(-0.1*b3*t))"
st<-c(b1=1, b2=1, b3=1)
low<-c(-Inf,-Inf, -Inf)
up<-c(2, 2, 2)

cat("try nlxb\n")
anlxb2<-nlxb(hobmod, st, data=weeddata, lower=low, upper=up)
anlxb2

try(anls2p<-nls(hobmod, st, data=weeddata,  lower=low, upper=up, algorithm='port'))
summary(anls2p)
anls2p$m$deviance()

aLM2<-nlsLM(hobmod, st, data=weeddata, lower=low, upper=up)
summary(aLM2)
aLM2$m$deviance()