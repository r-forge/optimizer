## @knitr nlsdata.R
# try different ways of supplying data to R nls stuff
y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
50.156, 62.948, 75.995, 91.972)
tt <- seq_along(y) # for testing

mydata <- data.frame(y = y, tt = tt)

hobsc <- y ~ 100*b1/(1 + 10*b2 * exp(-0.1 * b3 * tt))

ste <- c(b1 = 2, b2 = 5, b3 = 3)

# let's try finding the variables

findmainenv <- function(formula, prm) {
   vn <- all.vars(formula)
   pnames <- names(prm)
   ppos <- match(pnames, vn)
   datvar <- vn[-ppos]
   cat("Data variables:")
   print(datvar)
   cat("Are the variables present in the current working environment?\n")
   for (i in seq_along(datvar)){
       cat(datvar[[i]]," : present=",exists(datvar[[i]]),"\n")
   }
}

findmainenv(hobsc, ste)
# ===============================


# let's try finding the variables in dotargs



finddotargs <- function(formula, prm, ...) {
   dots <- list(...)
   cat("dots:")
   print(dots)
   cat("names in dots:")
   dtn <- names(dots)
   print(dtn)
   vn <- all.vars(formula)
   pnames <- names(prm)
   
   ppos <- match(pnames, vn)
   datvar <- vn[-ppos]
   cat("Data variables:")
   print(datvar)
   cat("Are the variables present in the dot args?\n")
   for (i in seq_along(datvar)){
       dname <- datvar[[i]]
       cat(dname," : present=",(dname %in% dtn),"\n")
   }
}

finddotargs(hobsc, ste, y=y, tt=tt)
# ===============================

nlsquiet <- nls(formula=hobsc, start=ste)
print(nlsquiet)
#- OK
nlsdots <- nls(formula=hobsc, start=ste, y=y, tt=tt)
print(nlsdots)
#- OK
nlsframe <- nls(formula=hobsc, start=ste, data=mydata)
print(nlsframe)
#- OK

library(nlsr)
nlsrquiet <- nlxb(formula=hobsc, start=ste)
print(nlsrquiet)
#- OK
test <- try(nlsrdots <- nlxb(formula=hobsc, start=ste, y=y, tt=tt))
  if (class(test) != "try-error") { print(nlsrdots) } else {cat("Try error\n") }
#- Note -- does NOT work -- do we need to specify the present env. in nlfb for y, tt??
test2 <- try(nlsframe <- nls(formula=hobsc, start=ste, data=mydata))
if (class(test) != "try-error") {print(nlsframe) } else {cat("Try error\n") }
#- OK


library(minpack.lm)
nlsLMquiet <- nlsLM(formula=hobsc, start=ste)
print(nlsLMquiet)
#- OK
## Dotargs
##?? nlsLMdots <- nlsLM(formula=hobsc, start=ste, y=y, tt=tt)
##?? print(nlsLMdots)
#-  Note -- does NOT work
## dataframe
nlsLMframe <- nlsLM(formula=hobsc, start=ste, data=mydata)
print(nlsLMframe)
#- OK

detach("package:nlsr", unload=TRUE)
library(nlmrt)
##?? nlxbquiet <- nlxb(formula=hobsc, start=ste)
##?? print(nlxbquiet)
#- Note -- does NOT work
##?? nlxbdots <- nlxb(formula=hobsc, start=ste, y=y, tt=tt)
##?? print(nlxbdots)
#- Note -- does NOT work
## dataframe
nlxbframe <- nlxb(formula=hobsc, start=ste, data=mydata)
print(nlxbframe)
#- OK

