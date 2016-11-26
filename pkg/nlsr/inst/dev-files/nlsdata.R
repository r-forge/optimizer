## nlsdata.R
# try different ways of supplying data to R nls stuff
y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
50.156, 62.948, 75.995, 91.972)
tt <- seq_along(y) # for testing

mydata <- data.frame(y = y, tt = tt)

hobsc <- y ~ 100*b1/(1 + 10*b2 * exp(-0.1 * b3 * tt))

ste <- c(b1 = 2, b2 = 5, b3 = 3)

nlsquiet <- nls(formula=hobsc, start=ste)
print(nlsquiet)
#- OK

library(nlmrt)
nlxbquiet <- nlxb(formula=hobsc, start=ste)
print(nlxbquiet)
#- Note -- does NOT work

library(minpack.lm)
nlsLMquiet <- nlsLM(formula=hobsc, start=ste)
print(nlsLMquiet)
#- OK

## Dotargs
nlsdots <- nls(formula=hobsc, start=ste, y=y, tt=tt)
print(nlsdots)
#- OK

nlxbdots <- nlxb(formula=hobsc, start=ste, y=y, tt=tt)
print(nlxbdots)
#- Note -- does NOT work

nlsLMdots <- nlsLM(formula=hobsc, start=ste, y=y, tt=tt)
print(nlsLMdots)
#-  Note -- does NOT work

## dataframe
nlsframe <- nls(formula=hobsc, start=ste, data=mydata)
print(nlsframe)
#- OK

nlxbframe <- nlxb(formula=hobsc, start=ste, data=mydata)
print(nlxbframe)
#- OK

nlsLMframe <- nlsLM(formula=hobsc, start=ste, data=mydata)
print(nlsLMframe)
#- OK

