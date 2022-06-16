# Meyer_10 function as a least squares problem
y <- c(34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030,
         6005, 5147, 4427, 3820, 3307, 2872)
m <- 16
t <- 45 + 5 * (1:m)
df <- data.frame(t, y)
modl <- y ~ x1 * exp(x2/(t + x3))
modlp <- y ~ exp(x2/(t + x3))
library(minpack.lm)
library(nlsr)
cat("nls:\n")
anls<-nls(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anls) # fails singular gradient
cat("NOTE: modlp NOT modl\n")
anlsp<-nls(formula=modlp, start=c(x2=1, x3=1), data=df, algorithm="plinear", trace=TRUE)
summary(anlsp) # fails singular gradient
anlxb<-nlxb(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlxb) # gets to min but slowly
anlxb
anlsLM<-nlsLM(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsLM) ## not near min
library(nlsj)
anlsj<-nlsj(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsj) # NOT good
library(nlsralt)
anlxbx<-nlxbx(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlxbx) # gets to min but slowly (using unscaled)
anlxbx


