options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test chen-x.f ...\n")

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

chen.g <- function(x) {
   res <- chen.res(x)
   jj<-chen.jac(x)
   gg<- 2.0 * jj %*% as.vector(res)
   #return(gg)
}

chen.res <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
}

chen.jac <- function(x) {
   n<-length(x)
   v <- log(x) + exp(x)
   jj<-matrix(0.0, n, n)
   for (i in 1:n) {
     jj[i,i] <- 0.5 * (1.0/x[i] + exp(x[i])) * (1.0-v[i]/sqrt(v[i]^2 + 5e-04))
   } 
   jj #return(jj)
}

p0 <- rexp(50)

cat("check fn at p0:")
fval<-chen.f(p0)
print(fval)
cat("check gr at p0:")
gval<-chen.g(p0)
print(gval)

# NOTE: lower bound to prevent non-computable function

system.time(anschenf <- optimx(par=p0, fn=chen.f, lower=0.001, 
   control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

optansout(anschenf, filename="./anschenf.txt")

system.time(anschenG <- optimx(par=p0, fn=chen.f, gr=chen.g, lower=0.001, 
   control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

optansout(anschenG,filename="./anschenG.txt")


ans.optx<-transform(anschenf, pflag="F")
ans.optxG<-transform(anschenG, pflag="FG")

achen<-rbind(anschenf, anschenG)

idx<-order(as.vector(achen$method, mode="character"))
achen<-achen[idx, ]
achen
 

cat("====================== end chen_test ========================\n")


