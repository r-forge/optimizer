
options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

ntpar<-10 # set the size of the problem

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test broydt-x.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

broydt.g <- function(x) {
   n <- length(x)
   gg<-rep(NA,n) # gradient set to NA to start with
   gg[1]<- -2 + 2*x[1] + 4*x[3] + (6 - 2.0*x[1])*(1 - 2*x[2] + x[1]*(3 - 0.5*x[1])) -
                2*x[2]*(3 - 0.5*x[2])
   gg[2] <- -6 + 4*x[4] + 10*x[2] + 
                (6 - 2.0*x[2])*(1 - x[1] - 2*x[3] + x[2]*(3 - 0.5*x[2])) - 
                4*x[1]*(3 - 0.5*x[1]) - 2*x[3]*(3 - 0.5*x[3])
   tnm2 <- 3:(n-2)
   gg[tnm2]<- -6 + 4*x[tnm2-2] + 4*x[tnm2+2] + 
               10*x[tnm2] + (6 - 2.0*x[tnm2])*(1 - x[tnm2-1] - 2*x[tnm2+1] +
               x[tnm2]*(3 - 0.5*x[tnm2])) - 4*x[tnm2-1]*(3 - 0.5*x[tnm2-1]) - 
               2*x[tnm2+1]*(3 - 0.5*x[tnm2+1])
   gg[n-1] <- -6 + 4*x[n-3] + 10*x[n-1] + 
             (6 - 2.0*x[n-1])*(1 - x[n-2] - 2*x[n] + x[n-1]*(3 - 0.5*x[n-1])) - 
             4*x[n-2]*(3 - 0.5*x[n-2]) - 2*x[n]*(3 - 0.5*x[n])

   gg[n] <- -4 + 4*x[n-2] + 8*x[n] + 
             (6 - 2.0*x[n])*(1 - x[n-1] + x[n]*(3 - 0.5*x[n])) - 
             4*x[n-1]*(3 - 0.5*x[n-1])
   gg
}

p0 <- rnorm(ntpar, sd=1) # starting values n=ntpar
cat("Test: broydt, no gradient\n")
system.time(ans.broydtf <- optimx(par=p0, fn=broydt.f,control=list(maxit=25000,save.failures=TRUE,all.methods=TRUE)))[1]
optansout(ans.broydtf,filename="./ansbroydtf.txt")

cat("Test: broydt, with analytic gradient\n")
system.time(ans.broydtfg <- optimx(par=p0, fn=broydt.f, gr=broydt.g, 
             control=list(maxit=25000,save.failures=TRUE,all.methods=TRUE)))[1]
optansout(ans.broydtfg,filename="./ansbroydtfg.txt")

cat("================= end broydt_test ==================\n")
