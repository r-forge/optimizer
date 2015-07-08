## ox - interactive example

require(graphics)
require(optimz)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

####
##JN Since this is a single method, the details could have wrong structure
ans2<-optimx(c(-1.2,1), fr, grr, method = "BFGS", control=list(trace=2))
ans2
tmpin<-readline("cont?")
str(ans2)
tmpin<-readline("cont?")
det2<-attr(ans2,"details")
##JN cat("dim(det2)\n")
print(dim(det2))
tmpin<-readline("cont?")
coef(ans2)
print(dim(coef(ans2)))
tmpin<-readline("cont?")
ans2d<-as.data.frame(ans2)
print(ans2d)
attributes(ans2d)
tmpin<-readline("cont?")
ans2parval<-ans2[,1:(attr(ans2,"npar")+1)]
print(ans2parval)
class(ans2parval)

tmpin<-readline("cont?")
summary(ans2)
dim(summary(ans2))
tmpin<-readline("cont?")

