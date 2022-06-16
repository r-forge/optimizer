# Meyer_10 function as a least squares problem
y <- c(34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030,
         6005, 5147, 4427, 3820, 3307, 2872)
t <- 45 + 5 * (1:m)
df <- data.frame(t, y)
m <- 16
modl <- y ~ x1 * exp(x2/(t + x3))
library(minpack.lm)
library(nlsr)
cat("nls:\n")
anls<-nls(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anls)
anlxb<-nlxb(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlxb)
anlxb
anlsLM<-nlsLM(formula=modl, start=c(x1=1, x2=1, x3=1), data=df, trace=TRUE)
summary(anlsLM)
