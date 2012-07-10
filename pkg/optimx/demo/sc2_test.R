options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("optimx test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

neg.sc2.f <- function(x){
n <- length(x)
vec <- 1:n
-sum(vec * (exp(x) - x)) / 10
}

neg.sc2.g <- function(x){
n <- length(x)
vec <- 1:n
-vec * (exp(x) - 1) / 10
}

ntpar<-25

p0 <- runif(ntpar,min=-1, max=1)
system.time(anssc2f <- optimx(par=p0, fn=sc2.f, 
   control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

optansout(anssc2f,filename="./anssc2f.txt")
anssc2f


system.time(anssc2fneg <- optimx(par=p0, fn=neg.sc2.f, 
   control=list(maxit=2500, maximize=TRUE,save.failures=TRUE, all.methods=TRUE)))[1]

optansout(anssc2fneg,filename="./anssc2negf.txt")
anssc2fneg


system.time(anssc2g <- optimx(par=p0, fn=sc2.f, gr=sc2.g,
   control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

optansout(anssc2g,filename="./anssc2g.txt")
anssc2g

system.time(anssc2gneg <- optimx(par=p0, fn=neg.sc2.f, gr=neg.sc2.g,
   control=list(maxit=2500, maximize=TRUE,save.failures=TRUE,all.methods=TRUE)))[1]

optansout(anssc2gneg,filename="./anssc2gneg.txt")
anssc2gneg

cat("================== end sc2_test ===================\n")

